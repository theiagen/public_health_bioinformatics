/*
 * Makes documentation tables resizable by dragging a column border, while
 * keeping them fully responsive to the CURRENT viewport width.
 *
 * Sizing strategy:
 *   - By default a table just uses CSS: `width: 100%` + `table-layout: auto`
 *     with a `min-width` on every cell (see extra.css). The browser therefore
 *     fills the current viewport live and only shows a horizontal scrollbar once
 *     the columns would have to shrink below that minimum — no JavaScript, no
 *     stale "width at load time", and it works inside inactive tabs too.
 *   - When the user drags a column border, that table switches to fixed layout
 *     with explicit widths so the drag is precise (neighbour-absorbing: widening
 *     a column narrows the one to its right). On release the widths are stored as
 *     PERCENTAGES, so the manually-sized table stays responsive to later
 *     viewport changes.
 *
 * It cooperates with the other table scripts:
 *   - it never adds a `class` to the <table> (it uses a data- attribute), so
 *     tablesort.js still enhances it and the CSS selectors still match,
 *   - handle drags call stopPropagation so a resize never triggers a sort,
 *   - the search box in table-search.js keeps working (it finds the table via a
 *     recursive querySelector, regardless of the new wrapper).
 */

const MIN_COLUMN_WIDTH = 70; // px — manual-drag floor (matches [data-manual] min-width in CSS)

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

  // a handle sits on each internal column border (not the table's outer edge)
  for (let i = 0; i < headers.length - 1; i++) {
    makeColumnResizable(headers[i], headers[i + 1], table, headers);
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
