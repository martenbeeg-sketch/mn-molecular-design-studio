# Institute teaching access — three students and one instructor

## Open in a browser

**Teaching URL:** <https://pc002975.marionegri.it:18502/>

The equivalent local-IP URL is <https://172.18.24.170:18502/>.
Use `https://`, not `http://`. Streamlit and the calculations run on PC002975;
the other computers only need a browser. This is the isolated, supervised
teaching pilot, not the existing research production deployment.

Access was verified from PC004032 and PC004033 using the `compute` operator
account. Four independent Streamlit sign-in sessions rendered concurrently on
each client, over certificate-verified TLS WebSockets. Four simultaneous
authenticated account/API reads were also checked separately. These checks are
not a completed authenticated-browser lesson or a GPU workload acceptance test.

## One-time certificate setup

No institute-issued server certificate was available in this setup. The pilot
uses a separate, short-lived teaching certificate. Before entering passwords,
have the instructor or institute IT approve and import its **public CA
certificate** into each participating computer's browser/user trust store:

`/home/user/Downloads/MDS-Teaching-Pilot-CA.crt`

Verify this SHA-256 certificate fingerprint through the instructor:

```text
D4:35:14:FB:78:98:84:A2:0E:E9:98:D0:EE:93:B7:CD:06:D7:3E:2D:7D:9A:EE:C4:0C:38:4D:1D:5B:FA:2C:1F
```

Only distribute this public `.crt` file, never any `.key` file or the private
teaching-state directory. Do not bypass certificate warnings, disable TLS
verification, or install a CA contrary to institute policy. An institute-issued
certificate for the hostname can replace this certificate without changing the
application or its data. No workstation's system/browser trust store was
automatically changed during deployment.

The CA is name-constrained to this hostname and IP. The server certificate
expires **2026-12-11 20:58:11 UTC**; the CA expires 2027-02-09. Renewal is manual,
not automatic. Do not schedule a lesson past the server expiry without renewal.

## Accounts and first lesson

- Instructor: `teaching_instructor`, application admin.
- Students: `teaching_student_one`, `teaching_student_two`,
  `teaching_student_three`, ordinary application users.
- The instructor's password and three individual temporary student passwords
  are in the private operator file:
  `/home/user/programs/deployments/mn-molecular-design-studio/teaching-pilot-current/teaching-state/credentials.json`.
- Give each student only their assigned temporary password through an approved
  private channel, not the complete credentials file. All three accounts require
  a password change at first sign-in; earlier automated student sessions were
  revoked. Their chosen new passwords are not written back to this operator file.
- These are four distinct teaching seats, not verified real-person identities.
  Record which student receives each seat; do not share the instructor account
  or the workstation's `compute`/`user` login or Docker access.

For the first supervised lesson, use the disposable shared teaching project:
browse the supplied molecules/public structure, run the prescribed native CPU
sequence-liability assessment, inspect its provenance, review and export results.
Do not upload important research data or start GPU/design campaigns. The pilot
has **four CPU cores and 8 GiB RAM total**, not per student. Persistent disk
quotas and per-student fair scheduling are not implemented. HTTP upload and
Streamlit message limits are 10 MB; these are not cumulative storage quotas.

## Network and operator boundaries

- HTTPS binds only `172.18.24.170:18502`, not all interfaces. The permitted client
  subnet is the directly connected `172.18.24.0/22`, which includes both tested
  compute PCs. It is not assumed to cover every institute VLAN or VPN.
- Other source subnets receive HTTP 403, even with a forged `X-Forwarded-For`.
  Obtain the actual approved subnet if a student's institute machine is outside
  this range; do not replace the allowlist with all private networks or `0/0`.
- The API stays loopback-only at `127.0.0.1:18766`; `/api/*` is not forwarded by
  the LAN entry point. Local instructor access remains
  <http://127.0.0.1:18502/> on PC002975 only.
- Caddy 2.11.9 terminates TLS and proxies only to the existing loopback teaching
  UI. Its admin API, automatic certificate requests, HTTP listener/redirects,
  HTTP/3 and automatic trust-store installation are disabled. No firewall or
  institute routing policy was modified.
- The gateway runs as the ordinary user in `mds-teaching-lan.service`, with
  256 MiB memory and one CPU core limits. Its config, TLS keys, pinned binary and
  checksum receipt live under `lan-gateway/` in the current private teaching
  release directory. The download was verified against the official release's
  SHA-512 checksum; this is not an independent release-signature verification.
- Caddy archive SHA-256:
  `5b8cb7faf8726a2028ad96a66e54e6e81e8f361bc6175e5d22bbd997cc5177fb`.
- `scripts/teaching_lan_gateway.py` renders the constrained configuration and
  rejects wildcard/public listeners, mismatched subnets and shared-readable
  TLS private keys. Validate the result with the pinned Caddy before activation.
- No service is enabled for automatic reboot startup. The host and pilot must
  remain running for the class. The private teaching gateway recovery copy
  includes its certificate/key/config/unit; do not put that copy in Git or email.

Stop LAN access without touching research services:

```bash
systemctl --user stop mds-teaching-lan.service
```

To start the existing retained pilot after an operator stop, start
`mds-teaching-pilot-20261111`, verify the teaching loopback relay is running, then
start `mds-teaching-lan.service`. A full workstation reboot also requires the
relay to be started from its recorded command; do not rerun account bootstrap
or substitute production state. The gateway alone cannot restart the pilot.

## Evidence and remaining limits

`TEACHING_LAN_EVIDENCE.tar.gz` contains sanitized network/account receipts, test
results, verification scripts and the source-denial/TLS checks. It contains no
passwords, TLS private keys, CA keys, databases, sessions or full recovery copies.
Its SHA-256 is
`f053ee004b4985f723ca3a3b165fbf1f832bc51392e6cc67fedf8f1beaf2b39b`.

The focused isolated suite passed **73 tests**. Network checks passed from all
three PCs. The original analysis run `20261111-190956-7720bdfa` remains the only
analysis job. The pilot's account/data restore was checked again after adding
the third student: 63 files/links and four SQLite databases passed offline
verification, including all four accounts and their first-login requirements.
The private gateway files have a separate recovery copy.

Research production services, authentication and data remain unchanged. The
production rollout, GPU workflows, cleanup decisions, historical-job
reconciliation and broader acceptance gates in `STUDENT_READINESS.md` remain
open. Successful class access is not a scientific or production certification.
