export default function ({data, parentElement, setTriggerValue}) {
    const workspace = parentElement.closest(`.st-key-${data.workspace}`);
    const row = workspace?.querySelector('[data-testid="stHorizontalBlock"]');
    const columns = row && [...row.children].filter(child => child.dataset.testid === 'stColumn');
    if (columns?.length !== data.widths.length) return;

    const focused = [...row.querySelectorAll(':scope > .studio-pane-divider')].indexOf(document.activeElement);
    row.studioPaneDividerCleanup?.();
    row.classList.add('studio-pane-split');
    const weights = [...data.widths];
    // Each divider changes the pane on its side of the flexible column.
    const dividerPanes = Object.keys(data.controls).map(Number);
    let drag;
    let pending = false;
    const handles = dividerPanes.map(pane => {
        const {label, limits} = data.controls[pane];
        const handle = document.createElement('div');
        handle.className = 'studio-pane-divider';
        handle.tabIndex = 0;
        handle.title = `Resize ${label.toLowerCase()}`;
        handle.setAttribute('role', 'separator');
        handle.setAttribute('aria-label', handle.title);
        handle.setAttribute('aria-orientation', 'vertical');
        handle.setAttribute('aria-valuemin', limits[0]);
        handle.setAttribute('aria-valuemax', limits[1]);
        handle.setAttribute('aria-valuenow', data.widths[pane]);
        row.appendChild(handle);
        return handle;
    });

    function positionHandles() {
        const panes = columns.map(column => column.getBoundingClientRect());
        const left = row.getBoundingClientRect().left;
        handles.forEach((handle, index) => {
            const first = panes[index];
            const second = panes[index + 1];
            handle.hidden = Math.abs(first.top - second.top) > 2;
            handle.style.left = `${(first.right + second.left) / 2 - left}px`;
        });
    }

    function preview(index, delta) {
        const pane = dividerPanes[index];
        const [minimum, maximum] = data.controls[pane].limits;
        const value = Math.max(minimum, Math.min(maximum, data.widths[pane] + delta));
        weights[pane] = value;
        weights[data.flexible] = data.widths[data.flexible] + data.widths[pane] - value;
        const gap = getComputedStyle(row).columnGap;
        columns.forEach((column, pane) => {
            // Streamlit's native columns distribute the unused gap via flex-grow.
            column.style.flexBasis = `calc(${weights[pane]}% - ${gap})`;
        });
        handles[index].setAttribute('aria-valuenow', value);
        positionHandles();
    }

    function clearPreview() {
        columns.forEach(column => column.style.removeProperty('flex-basis'));
    }

    function commit(pane) {
        if (weights[pane] !== data.widths[pane]) {
            pending = true;
            setTriggerValue('resize', [pane, weights[pane]]);
        } else {
            clearPreview();
            positionHandles();
        }
    }

    handles.forEach((handle, index) => {
        const pane = dividerPanes[index];
        const direction = pane < data.flexible ? 1 : -1;
        handle.onpointerdown = event => {
            if (event.button !== 0 || drag || pending) return;
            event.preventDefault();
            handle.focus({preventScroll: true});
            handle.setPointerCapture(event.pointerId);
            drag = {handle, x: event.clientX, width: row.clientWidth};
            row.classList.add('studio-pane-dragging');
        };
        handle.onpointermove = event => {
            if (drag?.handle === handle) {
                preview(index, Math.round(direction * (event.clientX - drag.x) * 100 / drag.width));
            }
        };
        handle.onpointerup = event => {
            if (drag?.handle !== handle) return;
            handle.onpointermove(event);
            drag = undefined;
            row.classList.remove('studio-pane-dragging');
            commit(pane);
        };
        handle.onlostpointercapture = () => {
            if (drag?.handle !== handle) return;
            drag = undefined;
            row.classList.remove('studio-pane-dragging');
            weights[pane] = data.widths[pane];
            weights[data.flexible] = data.widths[data.flexible];
            handle.setAttribute('aria-valuenow', data.widths[pane]);
            clearPreview();
            positionHandles();
        };
        handle.onkeydown = event => {
            if (pending || (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight')) return;
            event.preventDefault();
            preview(index, direction * (event.key === 'ArrowRight' ? 1 : -1));
            commit(pane);
        };
    });

    const observer = new ResizeObserver(positionHandles);
    observer.observe(row);
    columns.forEach(column => observer.observe(column));
    positionHandles();
    handles[focused]?.focus({preventScroll: true});
    const cleanup = () => {
        observer.disconnect();
        handles.forEach(handle => handle.remove());
        row.classList.remove('studio-pane-split', 'studio-pane-dragging');
        clearPreview();
    };
    // A component update can re-run this initializer on the same Streamlit row.
    row.studioPaneDividerCleanup = cleanup;
    return () => {
        if (row.studioPaneDividerCleanup === cleanup) {
            cleanup();
            delete row.studioPaneDividerCleanup;
        }
    };
}
