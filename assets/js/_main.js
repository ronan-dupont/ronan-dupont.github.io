/* Retain responsive video embeds and the existing image gallery. */
$(function () {
  $('#main').fitVids();
  var reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  $("a[href$='.jpg'],a[href$='.jpeg'],a[href$='.JPG'],a[href$='.png'],a[href$='.gif']")
    .addClass('image-popup');
  $('.image-popup').magnificPopup({
    type: 'image',
    tLoading: 'Loading image #%curr%...',
    gallery: { enabled: true, navigateByImgClick: true, preload: [0, 1] },
    image: { tError: '<a href="%url%">Image #%curr%</a> could not be loaded.' },
    removalDelay: reducedMotion ? 0 : 160,
    mainClass: 'mfp-zoom-in',
    closeOnContentClick: true,
    midClick: true
  });
});
