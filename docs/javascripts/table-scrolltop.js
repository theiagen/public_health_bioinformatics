/*
 * Adds a floating "scroll to top of table" button to every content table that
 * scrolls internally (tall tables get their own vertical scrollbar). The button
 * fades in only once the table has been scrolled down past a threshold, and
 * clicking it eases the table's scroll position back to the top.
 *
 * Depends on window.mdTables (table-utils.js) for the shared page-load hook, the
 * content-table selector, and the `.md-table-wrap` scroll container — the same
 * wrapper table-resize.js uses, so the two features share one container.
 */
(function () {
  const { onPageLoad, getContentTables, getTableWrap } = window.mdTables;

  const SCROLLTOP_THRESHOLD = 200; // px scrolled before the button appears
  const SCROLLTOP_DURATION = 700; // ms — fixed scroll-to-top time regardless of distance

  // Up-chevron icon for the per-table "scroll to top" button.
  const SCROLLTOP_ICON =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
    '<path d="M12 8l-6 6 1.41 1.41L12 10.83l4.59 4.58L18 14z"/>' +
    '</svg>';

  // Animate an element's scrollTop back to 0 over a FIXED duration (so the trip
  // takes the same time from anywhere in the table).
  function smoothScrollToTop(el) {
    const start = el.scrollTop;
    if (start <= 0) {
      return;
    }

    // quintic ease-out: quick push-off, then decelerate the whole way for a very
    // soft landing (the softest part of the motion)
    function ease(t) {
      return 1 - Math.pow(1 - t, 5);
    }

    let startTime = null;
    function step(now) {
      if (startTime === null) {
        startTime = now;
      }
      const t = Math.min((now - startTime) / SCROLLTOP_DURATION, 1);
      el.scrollTop = start * (1 - ease(t));
      if (t < 1) {
        requestAnimationFrame(step);
      }
    }

    requestAnimationFrame(step);
  }

  // Add a floating "scroll to top of table" button. It lives in a non-scrolling
  // positioning layer wrapped around the scroll container, so it stays pinned to
  // the visible bottom-right corner while the table scrolls internally. It fades
  // in only once the table has been scrolled down past SCROLLTOP_THRESHOLD.
  function addScrollTopButton(wrap) {
    // Re-runs on instant navigation; skip if the layer is already in place.
    const parent = wrap.parentElement;
    if (parent && parent.classList.contains("md-table-scroll-region")) {
      return;
    }

    const region = document.createElement("div");
    region.className = "md-table-scroll-region";
    wrap.parentNode.insertBefore(region, wrap);
    region.appendChild(wrap);

    const button = document.createElement("button");
    button.type = "button";
    button.className = "table-scrolltop is-hidden";
    button.setAttribute("aria-label", "Scroll to top of table");
    button.innerHTML = SCROLLTOP_ICON;

    button.addEventListener("click", function () {
      smoothScrollToTop(wrap);
    });

    function updateVisibility() {
      button.classList.toggle("is-hidden", wrap.scrollTop < SCROLLTOP_THRESHOLD);
    }

    wrap.addEventListener("scroll", updateVisibility, { passive: true });
    updateVisibility();

    region.appendChild(button);
  }

  function addAllScrollTopButtons() {
    getContentTables().forEach(function (table) {
      addScrollTopButton(getTableWrap(table));
    });
  }

  onPageLoad(addAllScrollTopButtons);
})();
