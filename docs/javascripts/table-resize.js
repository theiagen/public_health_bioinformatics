/*
 * Makes every content table's columns resizable:
 *   - drag an internal header border to trade width between it and its right
 *     neighbour (the table's total width stays fixed), and
 *   - drag a handle on the table's far-right header border to widen the last
 *     column / the whole table (adding a horizontal scrollbar) or shrink it
 */
(function () {
  const { // these parameters and functions are set in table-utils.js
    MIN_COLUMN_WIDTH,
    onPageLoad,
    getHeaderRow,
    getContentTables,
    getTableWrap,
  } = window.mdTables;

  /* ===== column resize helper functions ===== */

  // Freeze the current (CSS-computed) column widths as explicit pixels and switch
  // to fixed layout so a drag can adjust them precisely.
  function pinColumnPercentages(table, headers) {
    // get widths of columns and calculate total width
    const widths = headers.map(th => th.getBoundingClientRect().width);
    const total = widths.reduce((sum, w) => sum + w, 0);

    // set the layout to fixed so the drag can move them
    // also freeze current width of the whole table
    table.style.tableLayout = "fixed";
    table.style.width = total + "px";

    // tell css to set each column's width to a percentage of the table
    headers.forEach((th, i) => {
      th.style.width = (widths[i] / total) * 100 + "%";
    });

    // add a [dataset-manual] true flag
    table.dataset.manual = "true";
  }

  // Like pinColumnPercentages, but freezes column widths in PIXELS instead of
  // percentages; used by the last-column edge handle, which changes the table's
  // TOTAL width
  function pinColumnWidthsPx(table, headers) {
    const widths = headers.map(th => th.getBoundingClientRect().width);
    const total = widths.reduce((sum, w) => sum + w, 0);

    table.style.tableLayout = "fixed";
    table.style.width = total + "px";

    // set each columns width in pixels
    headers.forEach((th, i) => {
      th.style.width = widths[i] + "px";
    });

    table.dataset.manual = "true";
  }

  // resizable inner column functionality
  function makeColumnResizable(header, nextHeader, table, headers) {
    // make resize icon
    const handle = document.createElement("span");
    handle.className = "col-resize-handle";
    handle.setAttribute("aria-hidden", "true");
    header.appendChild(handle);

    // initialize shared variables
    let startX;
    let startLeftPct;
    let startRightPct;

    // behavior when the resize handle is clicked
    function onPointerDown(event) {
      // stop text selection or native drag behavior
      event.preventDefault();
      event.stopPropagation();

      // grab current column percentage widths
      pinColumnPercentages(table, headers);

      // get width of table and starting position of cursor
      const tableWidth = table.getBoundingClientRect().width;
      startX = event.clientX;

      // calculate starting column width percentages for both sides
      // of the pointer
      startLeftPct =
        parseFloat(header.style.width) ||
        (header.getBoundingClientRect().width / tableWidth) * 100;
      startRightPct =
        parseFloat(nextHeader.style.width) ||
        (nextHeader.getBoundingClientRect().width / tableWidth) * 100;

      // see lines 533-536 in extra.css for active resizing appearance
      document.body.classList.add("col-resizing");
      // now watch where the pointer goes and when it is released
      document.addEventListener("pointermove", onPointerMove);
      document.addEventListener("pointerup", onPointerUp);
    }

    // behavior when the resize handle is dragged
    function onPointerMove(event) {
      // get current width
      const tableWidth = table.getBoundingClientRect().width;

      // calculate the new percentage for the column based on where the cursor moved to
      const deltaPct =
        ((event.clientX - startX) / tableWidth) * 100;

      // calculate minimum allowed width in percentage
      const minPct =
        (MIN_COLUMN_WIDTH / tableWidth) * 100;

      // calculate the initial percentage volume of both columns
      const totalPct = startLeftPct + startRightPct;

      // change the column percentages based on the new position
      let leftPct = startLeftPct + deltaPct;
      let rightPct = startRightPct - deltaPct;

      // enforce minimum column width
      if (leftPct < minPct) {
        leftPct = minPct;
        rightPct = totalPct - leftPct;
      }
      if (rightPct < minPct) {
        rightPct = minPct;
        leftPct = totalPct - rightPct;
      }

      // set the columns to be the new widths by percentage
      header.style.width = leftPct + "%";
      nextHeader.style.width = rightPct + "%";
    }

    // behavior when the resize handle is released
    function onPointerUp() {
      // stop listening for movement or release
      document.removeEventListener("pointermove", onPointerMove);
      document.removeEventListener("pointerup", onPointerUp);
      // remove resizing appearance class
      document.body.classList.remove("col-resizing");

      // tell everyone the table was resized
      document.dispatchEvent(new Event("table-resized"));
    }

    // activate listening for when resize handle is clicked
    handle.addEventListener("pointerdown", onPointerDown);

    // listen for when the resize handle is pressed-and-released
    // this prevents sorting whenever you resize
    handle.addEventListener("click", function (event) {
      event.stopPropagation();
    });
  }

  // resizable last column (which can change the entire table width)
  function makeLastColumnResizable(header, table, headers) {
    // make resize icon (specific for the edge)
    const handle = document.createElement("span");
    handle.className = "col-resize-handle col-resize-handle--edge";
    handle.setAttribute("aria-hidden", "true");
    header.appendChild(handle);

    // initialize shared variables
    let startX;
    let startWidth;
    let startTableWidth;

    // behavior when resize handle is clicked
    function onPointerDown(event) {
      event.preventDefault();
      event.stopPropagation();

      // freeze existing column widths in pixels
      pinColumnWidthsPx(table, headers);

      // get starting position, column width, and table width
      startX = event.clientX;
      startWidth = header.getBoundingClientRect().width;
      startTableWidth = table.getBoundingClientRect().width;

      // add the column resizing appearance class and start listening for movement
      document.body.classList.add("col-resizing");
      document.addEventListener("pointermove", onPointerMove);
      document.addEventListener("pointerup", onPointerUp);
    }

    // behavior when resize handle is dragged
    function onPointerMove(event) {
      // grow or shrink the last column, down to its minimum
      let width = startWidth + (event.clientX - startX);
      if (width < MIN_COLUMN_WIDTH) {
        width = MIN_COLUMN_WIDTH;
      }

      // the whole table changes by the same amount
      header.style.width = width + "px";
      table.style.width = (startTableWidth + (width - startWidth)) + "px";
    }

    // behavior when resize handle is released
    function onPointerUp() {
      // stop listening, remove class, tell everybody about it
      document.removeEventListener("pointermove", onPointerMove);
      document.removeEventListener("pointerup", onPointerUp);
      document.body.classList.remove("col-resizing");
      document.dispatchEvent(new Event("table-resized"));
    }

    // start listening for resize handle click
    handle.addEventListener("pointerdown", onPointerDown);

    // when the resize handle is clicked-and-released, don't sort
    handle.addEventListener("click", function (event) {
      event.stopPropagation();
    });
  }

  // When the reset button is pushed, restore original appearance;
  function resetResize(table) {
    const headerRow = getHeaderRow(table);
    if (!headerRow) {
      return;
    }

    const headers = Array.from(headerRow.cells);

    // reset styles
    table.style.tableLayout = "";
    table.style.width = "";

    // keep only the base minimum width and let the CSS size it all
    table.style.minWidth = `${headers.length * MIN_COLUMN_WIDTH}px`;

    headers.forEach((th) => {
      th.style.width = "";
    });

    // remove manual tag
    delete table.dataset.manual;
  }

  // return true when a table has been resized
  function isResizeActive(table) {
    return table.dataset.manual === "true";
  }

  /* ===== table resizing lifecycle ===== */

  // set up table resizing on a single table
  function enableResizableTable(table) {
    // set applicable tables up for resizing (done only once)
    if (table.dataset.resizable === "true") {
      return;
    }
    const headerRow = getHeaderRow(table);
    if (!headerRow) {
      return;
    }
    table.dataset.resizable = "true";

    const headers = Array.from(headerRow.cells);

    // set minimum column widths in pixels
    table.style.minWidth = `${headers.length * MIN_COLUMN_WIDTH}px`;

    // enable the table to have a scrollbar instead of falling off the page
    getTableWrap(table);

    // make all internal columns resizable
    for (let i = 0; i < headers.length - 1; i++) {
      makeColumnResizable(headers[i], headers[i + 1], table, headers);
    }

    // make the last column resizable (unique due to table width modification)
    if (headers.length > 0) {
      makeLastColumnResizable(headers[headers.length - 1], table, headers);
    }
  }

  // add table resizing to all tables
  function enableAllTables() {
    getContentTables().forEach(enableResizableTable);
  }

  // reset table resizing
  function resetAllTables() {
    document
      .querySelectorAll("article table[data-manual]")
      .forEach(resetResize);
  }

  // enable the reset to only appear once the table is finished resizing
  let resizeTimer;
  window.addEventListener("resize", () => {
    clearTimeout(resizeTimer);

    resizeTimer = setTimeout(() => {
      resetAllTables();
    }, 150);
  });

  // make this feature's reset + active-check available for table-search.js's reset button
  window.mdTables.resetResize = resetResize;
  window.mdTables.isResizeActive = isResizeActive;

  // enable all the tables when a page is loaded
  onPageLoad(enableAllTables);
})();
