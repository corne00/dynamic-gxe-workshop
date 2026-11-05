// Minimal inputs for the interactive figures (no external libraries).
// Both return an element with a `value` property that fires "input" events,
// so they work as Observable `viewof` cells.
function control(label, input, show) {
  const box = document.createElement("div");
  box.className = "dgxe-control";
  const text = document.createElement("label");
  text.textContent = label;
  const out = document.createElement("output");
  box.append(text, input, out);
  const update = () => {
    box.value = input.type === "range" ? +input.value : input.value;
    out.textContent = show(box.value);
  };
  input.addEventListener("input", () => {
    update();
    box.dispatchEvent(new Event("input", {bubbles: true}));
  });
  update();
  return box;
}

// Slider over a list of values; the value of the control is the index.
function slider(label, values, start = 0, format = String) {
  const input = document.createElement("input");
  Object.assign(input, {type: "range", min: 0, max: values.length - 1,
                        step: 1, value: start});
  return control(label, input, (i) => format(values[i]));
}

// Drop-down list; the value of the control is the selected option.
function select(label, options, start = options[0]) {
  const input = document.createElement("select");
  for (const o of options) input.add(new Option(o, o));
  input.value = start;
  return control(label, input, () => "");
}
