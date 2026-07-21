/*
 * Continuous-spin icon for the "toggle" admonition.
 *
 * A collapsible toggle admonition is a <details class="toggle"> whose summary
 * shows an arrow icon (drawn by CSS as summary::before). With pure CSS that icon
 * can only sit at one of two fixed angles, so opening and closing animate in
 * OPPOSITE directions. To make it always spin the same way (clockwise), we keep
 * a running rotation total per box and add 180° on every open/close.
 *
 * JavaScript can't set styles on a ::before pseudo-element directly, so we write
 * the running angle into a CSS custom property (variable) on the <details>, and
 * the stylesheet reads it:
 *     summary::before { transform: rotate(var(--toggle-rotation)); }
 *
 * Uses mdTables.onPageLoad (table-utils.js) for its page-load registration; it
 * has no other dependency on the table scripts.
 */
(function () {
  const { onPageLoad } = window.mdTables;

  function spinToggleAdmonitions() {
    document.querySelectorAll("details.toggle").forEach(function (details) {
      // Don't bind twice if this runs again (e.g. an instant-navigation re-render).
      if (details.dataset.spinBound === "true") {
        return;
      }
      details.dataset.spinBound = "true";

      // Match the starting angle to the current state: upright (0°) when closed,
      // flipped (180°) when it loads already open — so an open-by-default box looks
      // the same as one you opened yourself. There's no transition until we add the
      // "toggle-animate" class below, so this starting angle applies instantly (no
      // spin on page load); only later open/close clicks animate.
      let rotation = details.open ? 180 : 0;
      details.style.setProperty("--toggle-rotation", rotation + "deg");
      void details.offsetWidth; // commit the starting angle before enabling animation
      details.classList.add("toggle-animate");

      details.addEventListener("toggle", function () {
        rotation += 180; // always the same direction -> continuous clockwise
        details.style.setProperty("--toggle-rotation", rotation + "deg");
      });
    });
  }

  // Run on first load and after every instant-navigation page change.
  onPageLoad(spinToggleAdmonitions);
})();
