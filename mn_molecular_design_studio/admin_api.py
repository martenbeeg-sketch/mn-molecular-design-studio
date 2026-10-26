"""Account and project grants through the same services as Users & access."""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, ConfigDict, Field, SecretStr


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid")


class NewUser(Contract):
    username: str
    display_name: str = ""
    temporary_password: SecretStr
    project_ids: list[str] = Field(default_factory=list)


class UserState(Contract):
    disabled: bool = Field(strict=True)


class ResetPassword(Contract):
    temporary_password: SecretStr


class ProjectGrants(Contract):
    project_ids: list[str]


class Settings(Contract):
    results_directory: str | None = None
    reference_directory: str | None = None
    api_base_url: str | None = None
    catalog_highlight_days: int | None = Field(default=None, ge=1, le=365, strict=True)


def install(app, current_user, admin_user):
    from mn_molecular_design_studio import auth, workspace
    from mn_molecular_design_studio.cluster import config
    router = APIRouter(prefix="/api/v1/admin", tags=["administration"])

    def local_authority():
        if config.is_worker():
            raise HTTPException(409, "Manage accounts on the coordinator.")

    def call(function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except PermissionError as exc:
            raise HTTPException(403, str(exc)) from exc
        except (OSError, ValueError, TypeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    def get_user(user_id):
        user = next((row for row in auth.list_users() if row["id"] == user_id), None)
        if not user:
            raise HTTPException(404, "User not found.")
        return {**user, "project_ids": auth.accessible_project_ids(user)}

    def checked_projects(values):
        result = []
        for value in dict.fromkeys(values):
            project = call(workspace.validate_project_id, value)
            if not project:
                raise HTTPException(422, "Choose an existing project.")
            result.append(project)
        return result

    @router.get("/users")
    def users(user=Depends(admin_user)):
        local_authority()
        return {"items": [{**row, "project_ids": auth.accessible_project_ids(row)} for row in auth.list_users()]}

    @router.post("/users", status_code=201)
    def create(payload: NewUser, user=Depends(admin_user)):
        local_authority()
        ids = checked_projects(payload.project_ids)
        created = call(auth.create_user, payload.username, payload.display_name,
                       payload.temporary_password.get_secret_value(), actor=user)
        for project in ids:
            call(auth.grant_project_access, project, created["id"], actor=user)
        return get_user(created["id"])

    @router.get("/users/{user_id}")
    def account(user_id: str, user=Depends(admin_user)):
        local_authority()
        return get_user(user_id)

    @router.patch("/users/{user_id}")
    def state(user_id: str, payload: UserState, user=Depends(admin_user)):
        local_authority()
        get_user(user_id)
        call(auth.set_user_disabled, user_id, payload.disabled, actor=user)
        return get_user(user_id)

    @router.post("/users/{user_id}/password")
    def password(user_id: str, payload: ResetPassword, user=Depends(admin_user)):
        local_authority()
        get_user(user_id)
        call(auth.reset_user_password, user_id, payload.temporary_password.get_secret_value(), actor=user)
        return {"status": "password_reset", "must_change_password": True, "sessions_revoked": True}

    @router.put("/users/{user_id}/projects")
    def grants(user_id: str, payload: ProjectGrants, user=Depends(admin_user)):
        local_authority()
        previous = get_user(user_id)
        if previous["role"] == "admin":
            raise HTTPException(422, "Admin accounts already have access to all projects.")
        ids = checked_projects(payload.project_ids)
        # Validate the entire replacement before writing any grants.
        for project in ids:
            call(auth.grant_project_access, project, user_id, actor=user)
        for project in set(previous["project_ids"]) - set(ids):
            call(auth.revoke_project_access, project, user_id, actor=user)
        return get_user(user_id)

    @router.get("/settings")
    def settings(user=Depends(admin_user)):
        return workspace.load_settings()

    @router.patch("/settings")
    def update_settings(payload: Settings, user=Depends(admin_user)):
        updates = payload.model_dump(exclude_unset=True)
        if any(value is None for value in updates.values()):
            raise HTTPException(422, "Settings cannot be null.")
        return call(workspace.save_settings, {**workspace.load_settings(), **updates})

    app.include_router(router)
