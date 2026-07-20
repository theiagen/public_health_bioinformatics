/*
 * Makes every content table's columns resizable by dragging:
 *   - drag an internal column border to trade width between it and its right
 *     neighbour (the table's total width stays fixed), and
 *   - drag a handle on the table's far-right edge to widen the last column / the
 *     whole table (adding a horizontal scrollbar) or shrink it back.
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
 * Depends on window.mdTables (table-utils.js) for MIN_COLUMN_WIDTH, the page-load
 * hook, the content-table selector, the header-row lookup, and the shared
 * `.md-table-wrap` scroll container. It publishes mdTables.resetResize(table) and
 * mdTables.isResizeActive(table) so table-search.js's reset button can drop a
 * table's manual widths without knowing how resizing works.
 *
 * Cooperation notes: this never adds a `class` to the <table> (it flags manual
 * sizing with data-manual), so tablesort and the CSS selectors still match; and
 * handle drags call stopPropagation so a resize never triggers a sort.
 */
(function () {
  const {
    MIN_COLUMN_WIDTH,
    onPageLoad,
    getHeaderRow,
    getContentTables,
    getTableWrap,
  } = window.mdTables;

  /* ===== Column resizing ===== */

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
  // again. Published on mdTables so table-search.js's reset button can call it.
  function resetResize(table) {
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
  }

  // Whether a table currently has manual (dragged) column widths.
  function isResizeActive(table) {
    return table.dataset.manual === "true";
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

    // ensure the table sits in the shared scroll container (a horizontal
    // scrollbar appears once a drag widens it past the container)
    getTableWrap(table);

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
    getContentTables().forEach(enhanceTable);
  }

  function resetAllManualTables() {
    document
      .querySelectorAll("article table[data-manual]")
      .forEach(resetResize);
  }

  // A viewport resize drops every manual sizing so the CSS re-fits the columns.
  let resizeTimer;
  window.addEventListener("resize", () => {
    clearTimeout(resizeTimer);

    resizeTimer = setTimeout(() => {
      resetAllManualTables();
    }, 150);
  });

  // Publish this feature's reset + active-check for table-search.js's reset button.
  window.mdTables.resetResize = resetResize;
  window.mdTables.isResizeActive = isResizeActive;

  onPageLoad(enhanceAllTables);
})();
