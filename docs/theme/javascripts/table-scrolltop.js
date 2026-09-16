/*
 * Adds a "scroll to top" button that appears when a table scrolls past a threshold
 */
(function () {
  // shared functions from table-utils.js
  const { onPageLoad, getContentTables, getTableWrap } = window.mdTables;

  // scrolltop functionality parameters
  const SCROLLTOP_THRESHOLD = 200; // px scrolled before the button appears
  const SCROLLTOP_DURATION = 700; // ms — fixed scroll-to-top time regardless of distance

  // icon for the "scroll to top" button.
  const SCROLLTOP_ICON =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
    '<path d="M12 8l-6 6 1.41 1.41L12 10.83l4.59 4.58L18 14z"/>' +
    '</svg>';

  // the animation for the scroll to the top over a FIXED duration
  // animation will always take the same amount of time
  function smoothScrollToTop(scrollbar) {
    // get curent position in the table
    const start = scrollbar.scrollTop;
    if (start <= 0) {
      return;
    }

    // how the animation works
    function ease(t) {
      // t: animation progress from start (0) to end (1)
      // should start fast and slow at the end
      return 1 - Math.pow(1 - t, 5);
    }

    let startTime = null;
    // schedules the animation until its done
    function step(now) {
      // start the animation timer
      if (startTime === null) {
        startTime = now;
      }
      // converts the elapsed time into the 0 to 1 progress fraction
      const t = Math.min((now - startTime) / SCROLLTOP_DURATION, 1);
      // sets the scroll position each frame
      scrollbar.scrollTop = start * (1 - ease(t));
      if (t < 1) {
        // keep scheduling the frame change until it finishes
        requestAnimationFrame(step);
      }
    }

    // activate the animation
    requestAnimationFrame(step);
  }

  // adds the "scroll to top of table" button. it appears only after the table
  // is scrolled past SCROLLTOP_THRESHOLD (in px)
  function addScrollTopButton(wrap) {
    // run upon navigation but skip if the layer is already in place.
    const parent = wrap.parentElement;
    if (parent && parent.classList.contains("md-table-scroll-region")) {
      return;
    }

    // add the scroll region class outside of the scroll area but relative to it
    const region = document.createElement("div");
    region.className = "md-table-scroll-region";
    wrap.parentNode.insertBefore(region, wrap);
    region.appendChild(wrap);

    // style the button and class features
    const button = document.createElement("button");
    button.type = "button";
    button.className = "table-scrolltop is-hidden";
    button.setAttribute("aria-label", "Scroll to top of table");
    button.innerHTML = SCROLLTOP_ICON;

    // listen for when the button is clicked to scroll to the top
    button.addEventListener("click", function () {
      smoothScrollToTop(wrap);
    });

    // hide it if you haven't scrolled past the threshold
    function updateVisibility() {
      button.classList.toggle("is-hidden", wrap.scrollTop < SCROLLTOP_THRESHOLD);
    }

    // listen for scrolling position
    wrap.addEventListener("scroll", updateVisibility, { passive: true });
    updateVisibility();

    // add the button
    region.appendChild(button);
  }

  // add the scroll to top button to all tables
  function addAllScrollTopButtons() {
    getContentTables().forEach(function (table) {
      addScrollTopButton(getTableWrap(table));
    });
  }

  // scroll the entire table to the top without using the button
  // primarily for table reset button
  function scrollTableToTop(table) {
    smoothScrollToTop(getTableWrap(table));
  }

  // when the page is loaded, add the scrollbars
  onPageLoad(addAllScrollTopButtons);

  // allow the table reset button to scroll the table to top
  window.mdTables.scrollTableToTop = scrollTableToTop;
})();
