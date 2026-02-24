// hpv_toggle.js

document.addEventListener("click", function (e) {
  const target = e.target;

  // ------------------------------------------------------------
  // HPV mapping pill animation (existing behavior)
  // ------------------------------------------------------------
  if (target.classList && target.classList.contains("hpv-pill")) {
    target.classList.add("hpv-pill-animate");
    setTimeout(function () {
      target.classList.remove("hpv-pill-animate");
    }, 220);

    // Prevent any other document click handlers from acting on this click
    // (useful if hpv-pill ever ends up inside other clickable containers)
    if (e.stopImmediatePropagation) e.stopImmediatePropagation();
    return;
  }

  // ------------------------------------------------------------
  // Open bucket when user clicks an input-preproc pill in the card
  // ------------------------------------------------------------
  const pillBtn = target.closest ? target.closest(".pp-card-pill-btn") : null;
  if (pillBtn) {
    const step = pillBtn.getAttribute("data-step");
    if (step && window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue("pp_open", step, { priority: "event" });
    }
    return;
  }

  // ------------------------------------------------------------
  // Preprocessing bucket pick buttons
  // ------------------------------------------------------------
  const pickBtn = target.closest ? target.closest(".pp-pick") : null;
  if (pickBtn) {
    if (pickBtn.disabled || pickBtn.classList.contains("is-disabled")) return;

    const kind = pickBtn.dataset.kind;
    const value = pickBtn.dataset.value;

    if (window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue(
        "pp_pick",
        { kind: kind, value: value, nonce: Date.now() },
        { priority: "event" }
      );
    }
    return;
  }

  // ------------------------------------------------------------
  // Bucket close button
  // ------------------------------------------------------------
  if (target.classList && target.classList.contains("pp-bucket-close")) {
    if (window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue("pp_bucket_close", Date.now(), { priority: "event" });
    }
    return;
  }
});

// Smooth scroll helper for bucket opening inside the wizard modal
document.addEventListener("DOMContentLoaded", function () {
  if (!window.Shiny) return;

  Shiny.addCustomMessageHandler("pp_scroll_to", function (msg) {
    const el = document.getElementById(msg.id);
    if (!el) return;
    try {
      el.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (e) {
      // no-op
    }
  });
});

document.addEventListener("keydown", function (e) {
  const pillBtn = e.target.closest ? e.target.closest(".pp-card-pill-btn") : null;
  if (!pillBtn) return;

  if (e.key === "Enter" || e.key === " ") {
    const step = pillBtn.getAttribute("data-step");
    if (step && window.Shiny && Shiny.setInputValue) {
      Shiny.setInputValue("pp_open", step, { priority: "event" });
    }
    e.preventDefault();
  }
});
