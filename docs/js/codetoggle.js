// Code visibility switch: off | R | Python (kept in localStorage).
(function () {
  const KEY = "dgxe-code";
  const get = () => { try { return localStorage.getItem(KEY) || "off"; } catch (e) { return "off"; } };
  const set = v => { try { localStorage.setItem(KEY, v); } catch (e) {} };
  const root = document.documentElement;
  function apply(mode) {
    root.classList.toggle("code-off", mode === "off");
    root.classList.toggle("code-on", mode !== "off");
    document.querySelectorAll(".code-switch button").forEach(b =>
      b.classList.toggle("active", b.dataset.mode === mode));
    if (mode !== "off") {
      // select the matching tab in every language-grouped tabset
      document.querySelectorAll('.panel-tabset[data-group="language"] > .nav-tabs .nav-link')
        .forEach(a => { if (a.textContent.trim() === mode && !a.classList.contains("active")) a.click(); });
    }
  }
  apply(get());                       // before first paint: no flash of code
  document.addEventListener("DOMContentLoaded", () => {
    document.querySelectorAll(".code-switch button").forEach(b =>
      b.addEventListener("click", () => { set(b.dataset.mode); apply(b.dataset.mode); }));
    apply(get());
    window.addEventListener("load", () => apply(get()));
  });
})();
