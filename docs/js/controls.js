// Minimal inputs for the interactive figures (no external libraries).
// Both return an element with a `value` property that fires "input" events,
// so they work as Observable `viewof` cells.
function control(label, input, read, show) {
  const box = document.createElement("div");
  box.className = "dgxe-control";
  const text = document.createElement("label");
  text.textContent = label;
  const out = document.createElement("output");
  box.append(text, input, out);
  const update = () => {
    box.value = read(input.value);
    out.textContent = show(box.value);
  };
  input.addEventListener("input", () => {
    update();
    box.dispatchEvent(new Event("input", {bubbles: true}));
  });
  update();
  return box;
}

function rangeInput(n, start) {
  const input = document.createElement("input");
  Object.assign(input, {type: "range", min: 0, max: n - 1, step: 1,
                        value: start});
  return input;
}

// Slider over a list of values; the value of the control is the index.
function slider(label, values, start = 0, format = String) {
  return control(label, rangeInput(values.length, start), (i) => +i,
                 (i) => format(values[i]));
}

// Slider over a list of values; the value of the control is the value.
function knob(label, values, start = 0, format = String) {
  return control(label, rangeInput(values.length, values.indexOf(start)),
                 (i) => values[i], format);
}

// Slider over indices with a play button that steps through them.
function player(label, values, start = 0, format = String, delay = 60) {
  const box = slider(label, values, start, format);
  const input = box.querySelector("input");
  const button = document.createElement("button");
  button.type = "button";
  button.className = "dgxe-play";
  button.textContent = "\u25B6";
  let timer = null;
  const stop = () => {
    clearInterval(timer);
    timer = null;
    button.textContent = "\u25B6";
  };
  button.addEventListener("click", () => {
    if (timer) return stop();
    if (+input.value >= values.length - 1) input.value = 0;
    button.textContent = "\u275A\u275A";
    timer = setInterval(() => {
      if (+input.value >= values.length - 1) return stop();
      input.value = +input.value + 1;
      input.dispatchEvent(new Event("input"));
    }, delay);
  });
  box.insertBefore(button, input);
  return box;
}

// Drop-down list; the value of the control is the selected option.
function select(label, options, start = options[0]) {
  const input = document.createElement("select");
  for (const o of options) input.add(new Option(o, o));
  input.value = start;
  return control(label, input, (v) => v, () => "");
}

// Remove the copy buttons of hidden code cells (outside the language tabs).
window.addEventListener("load", () => {
  for (const box of document.querySelectorAll(".code-copy-outer-scaffold")) {
    const pre = box.querySelector("pre");
    if (!box.closest(".tab-pane") && !(pre && pre.offsetHeight)) box.remove();
  }
});
