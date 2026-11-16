// Tell embedded pages (e.g. the LIMSport config builder) which light/dark
// scheme this site is showing, so they can match it.
//
// Why this is needed: an <iframe> shows a page from another website, and the
// browser does not let this page reach into it to change its styles. The two
// pages can only talk by sending messages (window.postMessage). The embedded
// page listens for {type: "limsport-theme", scheme: "light"|"slate"} and
// restyles itself.
//
// Opt an iframe in by adding the `data-sync-theme` attribute to it.
(function () {
  // The theme stamps the current scheme ("light" or "slate") on <body> as
  // data-md-color-scheme; check <html> too in case a theme version moves it.
  function currentScheme() {
    return (
      document.body.getAttribute("data-md-color-scheme") ||
      document.documentElement.getAttribute("data-md-color-scheme")
    );
  }

  function sendScheme(iframe) {
    const scheme = currentScheme();
    // contentWindow is null until the iframe exists in the page
    if (scheme && iframe.contentWindow) {
      // "*" = deliver whatever site the iframe is showing; the message holds
      // no private data, and the receiver only accepts the two known values
      iframe.contentWindow.postMessage({ type: "limsport-theme", scheme: scheme }, "*");
    }
  }

  function sendToAll() {
    document.querySelectorAll("iframe[data-sync-theme]").forEach(sendScheme);
  }

  // On every page load (including the theme's "instant" page switches, which
  // swap content without a full reload): send now, in case the iframe already
  // loaded, and again when it finishes loading, since a message sent before
  // then is lost (loading="lazy" delays it until scrolled near).
  window.mdTables.onPageLoad(function () {
    document.querySelectorAll("iframe[data-sync-theme]").forEach(function (iframe) {
      sendScheme(iframe);
      iframe.addEventListener("load", function () { sendScheme(iframe); });
    });
  });

  // Re-send whenever the reader flips the light/dark toggle. A
  // MutationObserver calls us back when the watched attribute changes.
  const observer = new MutationObserver(sendToAll);
  [document.documentElement, document.body].forEach(function (el) {
    observer.observe(el, { attributes: true, attributeFilter: ["data-md-color-scheme"] });
  });
})();
