export default function ({data, parentElement, setTriggerValue}) {
    const workspace = parentElement.closest('.st-key-datasets-workspace');
    const row = workspace?.querySelector('[data-testid="stHorizontalBlock"]');
    const columns = row && [...row.children].filter(child => child.dataset.testid === 'stColumn');
    if (columns?.length !== 3) return;

    const focused = [...row.querySelectorAll(':scope > .studio-dataset-divider')].indexOf(document.activeElement);
    row.studioDatasetDividerCleanup?.();
    row.classList.add('studio-dataset-split');
    const [folders, details] = data.widths;
    const weights = [folders, 100 - folders - details, details];
    let drag;
    let pending = false;
    const handles = ['Resize dataset folders', 'Resize member details'].map((label, index) => {
        const handle = document.createElement('div');
        handle.className = 'studio-dataset-divider';
        handle.tabIndex = 0;
        handle.title = label;
        handle.setAttribute('role', 'separator');
        handle.setAttribute('aria-label', label);
        handle.setAttribute('aria-orientation', 'vertical');
        handle.setAttribute('aria-valuemin', data.bounds[index][0]);
        handle.setAttribute('aria-valuemax', data.bounds[index][1]);
        handle.setAttribute('aria-valuenow', data.widths[index]);
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
        const [minimum, maximum] = data.bounds[index];
        const value = Math.max(minimum, Math.min(maximum, data.widths[index] + delta));
        weights[index === 0 ? 0 : 2] = value;
        weights[1] = 100 - weights[0] - weights[2];
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

    function commit() {
        if (weights[0] !== folders || weights[2] !== details) {
            pending = true;
            setTriggerValue('resize', [weights[0], weights[2]]);
        } else {
            clearPreview();
            positionHandles();
        }
    }

    handles.forEach((handle, index) => {
        // Moving the right divider right makes the details pane narrower.
        const direction = index === 0 ? 1 : -1;
        handle.onpointerdown = event => {
            if (event.button !== 0 || drag || pending) return;
            event.preventDefault();
            handle.focus({preventScroll: true});
            handle.setPointerCapture(event.pointerId);
            drag = {handle, x: event.clientX, width: row.clientWidth};
            row.classList.add('studio-dataset-dragging');
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
            row.classList.remove('studio-dataset-dragging');
            commit();
        };
        handle.onlostpointercapture = () => {
            if (drag?.handle !== handle) return;
            drag = undefined;
            row.classList.remove('studio-dataset-dragging');
            weights[0] = folders;
            weights[1] = 100 - folders - details;
            weights[2] = details;
            handle.setAttribute('aria-valuenow', data.widths[index]);
            clearPreview();
            positionHandles();
        };
        handle.onkeydown = event => {
            if (pending || (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight')) return;
            event.preventDefault();
            preview(index, direction * (event.key === 'ArrowRight' ? 1 : -1));
            commit();
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
        row.classList.remove('studio-dataset-split', 'studio-dataset-dragging');
        clearPreview();
    };
    // A component update can re-run this initializer on the same Streamlit row.
    row.studioDatasetDividerCleanup = cleanup;
    return () => {
        if (row.studioDatasetDividerCleanup === cleanup) {
            cleanup();
            delete row.studioDatasetDividerCleanup;
        }
    };
}
