/* Accessible, progressively enhanced navigation. No background polling. */
(function () {
  'use strict';
  function disclosure(options) {
    var button = document.querySelector(options.button);
    var panel = document.querySelector(options.panel);
    if (!button || !panel) { return; }
    var compact = window.matchMedia(options.breakpoint);
    var open = false;
    function render() {
      var expanded = compact.matches && open;
      panel.hidden = compact.matches && !open;
      panel.classList.toggle('is-open', expanded);
      button.setAttribute('aria-expanded', String(expanded));
    }
    function close(restoreFocus) {
      open = false; render();
      if (restoreFocus) { button.focus(); }
    }
    button.addEventListener('click', function () { open = !open; render(); });
    document.addEventListener('keydown', function (event) {
      if (event.key === 'Escape' && open) { close(true); }
    });
    document.addEventListener('click', function (event) {
      if (open && !panel.contains(event.target) && !button.contains(event.target)) { close(false); }
    });
    panel.addEventListener('click', function (event) {
      if (event.target.closest('a') && compact.matches) { close(false); }
    });
    var onResize = function () {
      var restoreFocus = compact.matches && panel.contains(document.activeElement);
      open = false; render();
      if (restoreFocus) { button.focus(); }
    };
    if (compact.addEventListener) { compact.addEventListener('change', onResize); }
    else { compact.addListener(onResize); }
    render();
  }
  function setup() {
    disclosure({ button: '.studio-menu-toggle', panel: '#navigation-panel', breakpoint: '(max-width: 959px)' });
    disclosure({ button: '.profile-toggle', panel: '#author-links', breakpoint: '(max-width: 959px)' });
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', setup); }
  else { setup(); }
}());
