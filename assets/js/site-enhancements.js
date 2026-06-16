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

    document.addEventListener("click", function (event) {
      if (!wrapper.contains(event.target) && wrapper.classList.contains("is-open")) {
        wrapper.classList.remove("is-open");
        button.classList.remove("open");
        button.setAttribute("aria-expanded", "false");
        links.style.display = "none";
      }
    });
  }

  function setupCardReveal() {
    if (!("IntersectionObserver" in window)) {
      return;
    }

    var cards = document.querySelectorAll(".archive__item, .research-card, .cv-panel");
    if (!cards.length) {
      return;
    }

    var observer = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) {
          entry.target.classList.add("is-visible");
          observer.unobserve(entry.target);
        }
      });
    }, { rootMargin: "0px 0px -8% 0px", threshold: 0.08 });

    cards.forEach(function (card) {
      card.classList.add("reveal-on-scroll");
      observer.observe(card);
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    setupAuthorLinks();
    setupCardReveal();
  });
}());
