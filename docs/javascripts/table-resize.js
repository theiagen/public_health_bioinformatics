/* ============================================================================
 * PRIMER for non-JavaScript developers (shared by the three table-*.js files)
 * ============================================================================
 *
 * These scripts run in the reader's browser and enhance the plain HTML tables
 * that Markdown/Zensical produce — adding search, sorting, column resizing, and
 * a scroll-to-top button. A few browser/theme concepts recur throughout:
 *
 *  - The DOM is the live, in-memory tree of the page. You find elements with a
 *    CSS selector and then read or change them:
 *      document.querySelector(sel)      first element matching a selector
 *      document.querySelectorAll(sel)   all matching elements (loop with forEach)
 *      el.textContent                   the element's text
 *      el.classList.add / toggle        add or flip a CSS class
 *      el.style.display = "none"        set an inline style
 *      el.dataset.manual                read/write a data-manual="..." attribute
 *      el.addEventListener(evt, fn)     run fn when evt happens (click, input,
 *                                       scroll, pointerdown, ...)
 *
 *  - document$ is a Material for MkDocs "observable" (an event stream). We call
 *    document$.subscribe(fn) to run fn once when the page loads AND again after
 *    every "instant navigation" — Material swaps page content without a full
 *    reload, so an ordinary page-load event would not fire on later pages.
 *
 *  - MutationObserver watches an element for changes (e.g. an attribute being
 *    set) and calls a function when they happen.
 *
 *  - requestAnimationFrame(fn) asks the browser to run fn just before the next
 *    screen repaint (when the browser redraws the page). It batches rapid work
 *    (like filtering while the reader types) into one update per frame instead
 *    of one per keystroke.
 *
 *  - Function declarations are "hoisted": a `function foo() {}` can be called
 *    from code written above it, so these files define helpers in a readable
 *    order rather than strictly-before-use order.
 *
 * The three files cooperate: table-sort.js adds click-to-sort; table-resize.js
 * (this file) wraps each table in a scrolling container and adds resizable
 * columns plus the scroll-to-top button; table-search.js adds the search box
 * and reset button.
 * ============================================================================ */

/*
 * Table enhancements applied to every content table:
 *   - resizable columns (drag a column border), plus a handle on the table's
 *     far-right edge that widens the last column / the whole table (adding a
 *     horizontal scrollbar), and
 *   - a floating "scroll to top of table" button for tall tables that scroll
 *     internally.
 *
 * Sizing strategy:
 *   - By default a table just uses CSS: `width: 100%` + `table-layout: auto`
 *     with a `min-width` on every cell (see extra.css). The browser therefore
 *     fills the current viewport (the visible browser area) live and only shows
 *     a horizontal scrollbar once the columns would have to shrink below that
 *     minimum — no JavaScript, no stale "width at load time", and it works
 *     inside inactive tabs too.
 *   - When the user drags a column border, that table switches to fixed layout
 *     with explicit widths so the drag is precise (neighbour-absorbing: widening
 *     a column narrows the one to its right). On release the widths are stored as
 *     PERCENTAGES, so the manually-sized table stays responsive (responsive =
 *     adapts to the window size) to later viewport changes.
 *
 * It cooperates with the other table scripts:
 *   - it never adds a `class` to the <table> (it uses a data- attribute), so
 *     tablesort.js still enhances it and the CSS selectors still match,
 *   - handle drags call stopPropagation so a resize never triggers a sort,
 *   - the search box in table-search.js keeps working (it finds the table via a
 *     recursive querySelector, regardless of the new wrapper).
 */

/* ===== Constants ===== */

const MIN_COLUMN_WIDTH = 70; // px — manual-drag floor (matches [data-manual] min-width in CSS)
const SCROLLTOP_THRESHOLD = 200; // px scrolled before the scroll-to-top button appears
const SCROLLTOP_DURATION = 700; // ms — fixed scroll-to-top time regardless of distance

// Up-chevron icon for the per-table "scroll to top" button.
const SCROLLTOP_ICON =
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
  '<path d="M12 8l-6 6 1.41 1.41L12 10.83l4.59 4.58L18 14z"/>' +
  '</svg>';

/* ===== Column resizing ===== */

function getHeaderRow(table) {
  return (table.tHead && table.tHead.rows[0]) || table.rows[0] || null;
}

// Freeze the current (CSS-computed) column widths as explicit pixels and switch
// to fixed layout so a drag can adjust them precisely. Re-runnable: each drag
// re-measures the live widths (which may currently be percentages).
function pinColumnPercentages(table, headers) {
  const widths = headers.map(th => th.getBoundingClientRect().width);
  const total = widths.reduce((sum, w) => sum + w, 0);

  table.style.tableLayout = "fixed";
  table.style.width = total + "px";   // freeze current width

  headers.forEach((th, i) => {
    th.style.width = (widths[i] / total) * 100 + "%";
  });

  table.dataset.manual = "true";
}

// Like pinColumnPercentages, but freezes column widths in PIXELS instead of
// percentages. Used by the last-column edge handle, which changes the table's
// TOTAL width — percentages (relative to that changing total) would fight it.
function pinColumnWidthsPx(table, headers) {
  const widths = headers.map(th => th.getBoundingClientRect().width);
  const total = widths.reduce((sum, w) => sum + w, 0);

  table.style.tableLayout = "fixed";
  table.style.width = total + "px";

  headers.forEach((th, i) => {
    th.style.width = widths[i] + "px";
  });

  table.dataset.manual = "true";
}

function makeColumnResizable(header, nextHeader, table, headers) {
  const handle = document.createElement("span");
  handle.className = "col-resize-handle";
  handle.setAttribute("aria-hidden", "true");
  header.appendChild(handle);

  let startX;
  let startLeftPct;
  let startRightPct;

  function onPointerMove(event) {
    const tableWidth = table.getBoundingClientRect().width;

    const deltaPct =
      ((event.clientX - startX) / tableWidth) * 100;

    const minPct =
      (MIN_COLUMN_WIDTH / tableWidth) * 100;

    const totalPct = startLeftPct + startRightPct;

    let leftPct = startLeftPct + deltaPct;
    let rightPct = startRightPct - deltaPct;

    if (leftPct < minPct) {
      leftPct = minPct;
      rightPct = totalPct - leftPct;
    }

    if (rightPct < minPct) {
      rightPct = minPct;
      leftPct = totalPct - rightPct;
    }

    header.style.width = leftPct + "%";
    nextHeader.style.width = rightPct + "%";
  }

  function onPointerUp() {
    document.removeEventListener("pointermove", onPointerMove);
    document.removeEventListener("pointerup", onPointerUp);
    document.body.classList.remove("col-resizing");

    document.dispatchEvent(new Event("table-resized"));
  }

  handle.addEventListener("pointerdown", function (event) {
    event.preventDefault();
    event.stopPropagation();

    pinColumnPercentages(table, headers);

    const tableWidth = table.getBoundingClientRect().width;

    startX = event.clientX;

    startLeftPct =
      parseFloat(header.style.width) ||
      (header.getBoundingClientRect().width / tableWidth) * 100;

    startRightPct =
      parseFloat(nextHeader.style.width) ||
      (nextHeader.getBoundingClientRect().width / tableWidth) * 100;

    document.body.classList.add("col-resizing");
    document.addEventListener("pointermove", onPointerMove);
    document.addEventListener("pointerup", onPointerUp);
  });

  handle.addEventListener("click", function (event) {
    event.stopPropagation();
  });
}

// A handle on the LAST column's right edge grows or shrinks that column — and
// with it the whole table. Drag right to widen the table past its container (a
// horizontal scrollbar appears); drag left to shrink the whole table, down until
// the last column reaches its minimum width. The internal handles above keep the
// table width fixed and trade width between neighbours; this one instead changes
// the total width.
function makeLastColumnResizable(header, table, headers) {
  const handle = document.createElement("span");
  handle.className = "col-resize-handle col-resize-handle--edge";
  handle.setAttribute("aria-hidden", "true");
  header.appendChild(handle);

  let startX;
  let startWidth;       // last column width at drag start (px)
  let startTableWidth;  // table width at drag start (px)

  function onPointerMove(event) {
    // grow/shrink the last column 1:1 with the pointer, down to its minimum
    let width = startWidth + (event.clientX - startX);
    if (width < MIN_COLUMN_WIDTH) {
      width = MIN_COLUMN_WIDTH;
    }

    // the whole table changes by the same amount, so it can end up either wider
    // OR narrower than the container
    header.style.width = width + "px";
    table.style.width = (startTableWidth + (width - startWidth)) + "px";
  }

  function onPointerUp() {
    document.removeEventListener("pointermove", onPointerMove);
    document.removeEventListener("pointerup", onPointerUp);
    document.body.classList.remove("col-resizing");

    document.dispatchEvent(new Event("table-resized"));
  }

  handle.addEventListener("pointerdown", function (event) {
    event.preventDefault();
    event.stopPropagation();

    pinColumnWidthsPx(table, headers);

    startX = event.clientX;
    startWidth = header.getBoundingClientRect().width;
    startTableWidth = table.getBoundingClientRect().width;

    document.body.classList.add("col-resizing");
    document.addEventListener("pointermove", onPointerMove);
    document.addEventListener("pointerup", onPointerUp);
  });

  handle.addEventListener("click", function (event) {
    event.stopPropagation();
  });
}

// Drop a table's manual widths and let the responsive CSS take over sizing
// again. Exposed globally so table-search.js's reset button can call it.
window.resetTableResize = function resetTableResize(table) {
  const headerRow = getHeaderRow(table);

  if (!headerRow) {
    return;
  }

  const headers = Array.from(headerRow.cells);

  table.style.tableLayout = "";
  table.style.width = "";

  // Keep only the intrinsic minimum width. CSS takes over sizing again.
  table.style.minWidth = `${headers.length * MIN_COLUMN_WIDTH}px`;

  headers.forEach((th) => {
    th.style.width = "";
  });

  delete table.dataset.manual;
};

/* ===== Scroll-to-top button ===== */

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

/* ===== Enhancement + lifecycle ===== */

function enhanceTable(table) {
  if (table.dataset.resizable === "true") {
    return;
  }
  const headerRow = getHeaderRow(table);
  if (!headerRow) {
    return;
  }
  table.dataset.resizable = "true";

  const headers = Array.from(headerRow.cells);

  // Keep the table from shrinking below usable column widths
  table.style.minWidth = `${headers.length * MIN_COLUMN_WIDTH}px`;

  // wrap the table in a scroll container (reuse one if present)
  let wrap = table.parentElement;
  if (!wrap || !wrap.classList.contains("md-table-wrap")) {
    wrap = document.createElement("div");
    wrap.className = "md-table-wrap";
    table.parentNode.insertBefore(wrap, table);
    wrap.appendChild(table);
  }

  // floating "scroll to top of table" button, pinned to the scroll container
  addScrollTopButton(wrap);

  // a handle sits on each internal column border (not the table's outer edge)
  for (let i = 0; i < headers.length - 1; i++) {
    makeColumnResizable(headers[i], headers[i + 1], table, headers);
  }

  // a handle on the far-right edge widens the last column (and the whole table)
  if (headers.length > 0) {
    makeLastColumnResizable(headers[headers.length - 1], table, headers);
  }
}

function enhanceAllTables() {
  // enhance every content table (same set tablesort targets); the richer visual
  // styling in the CSS is still limited to `.searchable-table` tables
  document
    .querySelectorAll("article table:not([class])")
    .forEach(enhanceTable);
}

function resetAllManualTables() {
  document
    .querySelectorAll("article table[data-manual]")
    .forEach(window.resetTableResize);
}

// A viewport resize drops every manual sizing so the CSS re-fits the columns.
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);

  resizeTimer = setTimeout(() => {
    resetAllManualTables();
  }, 150);
});

// Run on load and on every instant-navigation page change.
if (typeof document$ !== "undefined" && document$.subscribe) {
  document$.subscribe(enhanceAllTables);
} else {
  document.addEventListener("DOMContentLoaded", enhanceAllTables);
}
