/* Small progressive enhancements for the academic website. */
(function () {
  "use strict";

  function setupAuthorLinks() {
    var wrapper = document.querySelector(".author__urls-wrapper");
    if (!wrapper) {
      return;
    }

    var button = wrapper.querySelector("button");
    var links = wrapper.querySelector(".author__urls");
    if (!button || !links) {
      return;
    }

    button.setAttribute("type", "button");
    button.setAttribute("aria-controls", "author-links");
    button.setAttribute("aria-expanded", "false");
    links.setAttribute("id", "author-links");

    button.addEventListener("click", function () {
      window.setTimeout(function () {
        var isOpen = button.classList.contains("open");
        wrapper.classList.toggle("is-open", isOpen);
        button.setAttribute("aria-expanded", isOpen ? "true" : "false");
      }, 0);
    });

    wrapper.addEventListener("keydown", function (event) {
      if (event.key === "Escape" && wrapper.classList.contains("is-open")) {
        button.click(); button.focus();
      }
    });

    document.addEventListener("click", function (event) {
      if (!wrapper.contains(event.target) && wrapper.classList.contains("is-open")) {
        wrapper.classList.remove("is-open");
        button.classList.remove("open");
        button.setAttribute("aria-expanded", "false");
        links.style.display = "none";
      }
    });
  }

  function setupNavigation() {
    var nav = document.getElementById("site-nav");
    if (!nav) { return; }
    var toggle = nav.querySelector("button");
    var overflow = document.getElementById("navigation-overflow");
    if (!toggle || !overflow) { return; }
    function syncState() {
      toggle.setAttribute("aria-expanded", String(!overflow.classList.contains("hidden")));
    }
    new MutationObserver(syncState).observe(overflow, { attributes: true, attributeFilter: ["class"] });
    nav.addEventListener("keydown", function (event) {
      if (event.key === "Escape" && !overflow.classList.contains("hidden")) {
        toggle.click(); toggle.focus();
      }
    });
    syncState();
  }

  function setupFooterSize() {
    var footer = document.querySelector(".page__footer");
    if (footer && "ResizeObserver" in window) {
      new ResizeObserver(function () {
        document.body.style.marginBottom = footer.getBoundingClientRect().height + "px";
      }).observe(footer);
    }
  }


  document.addEventListener("DOMContentLoaded", function () {
    setupAuthorLinks();
    setupNavigation();
    setupFooterSize();
  });
}());
