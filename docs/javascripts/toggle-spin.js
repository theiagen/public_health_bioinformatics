/*
 * Continuous-360°-spin icon for the "toggle" admonition icon
 */
(function () {
  // shared functions from table-utils.js
  const { onPageLoad } = window.mdTables;

  // spin the toggle icon
  function spinToggleAdmonitions() {
    // grab all toggle admonitions
    document.querySelectorAll("details.toggle").forEach(function (details) {
      // don't attach this feature if it's already attached
      if (details.dataset.spinBound === "true") {
        return;
      }
      details.dataset.spinBound = "true";

      // set starting position upright (0°) when closed, flipped (180°) when open
      let rotation = details.open ? 180 : 0;
      details.style.setProperty("--toggle-rotation", rotation + "deg");
      void details.offsetWidth; // commit the starting angle before enabling animation
      details.classList.add("toggle-animate");

      // listen for when the toggle is clicked
      details.addEventListener("toggle", function () {
        // spin the toggle!
        rotation += 180;
        details.style.setProperty("--toggle-rotation", rotation + "deg");
      });
    });
  }

  // add spinning toggles to all pages when loaded
  onPageLoad(spinToggleAdmonitions);
})();
