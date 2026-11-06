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

// Row of buttons, one per option; the value of the control is the selected one.
function buttons(label, options, start = options[0]) {
  const box = document.createElement("div");
  box.className = "dgxe-control";
  const text = document.createElement("label");
  text.textContent = label;
  const group = document.createElement("div");
  group.className = "dgxe-buttons";
  const items = options.map((o) => {
    const b = document.createElement("button");
    b.type = "button";
    b.textContent = o;
    b.addEventListener("click", () => {
      box.value = o;
      items.forEach((x) => x.classList.toggle("active", x === b));
      box.dispatchEvent(new Event("input", {bubbles: true}));
    });
    group.append(b);
    return b;
  });
  box.append(text, group);
  box.value = start;
  items[options.indexOf(start)].classList.add("active");
  return box;
}

// Optimum and upper limit sliders for several genotypes, with a reset button.
// specs: [{name, colour, optimum, limit}]; the value is a copy of the specs.
function genotypeControls(specs) {
  const box = document.createElement("div");
  box.className = "dgxe-genotypes";
  const state = specs.map((s) => ({...s}));
  const rows = specs.map((spec, i) => {
    const row = document.createElement("div");
    row.className = "dgxe-genotype";
    const name = document.createElement("span");
    name.className = "dgxe-name";
    name.style.color = spec.colour;
    name.textContent = spec.name;
    const add = (label, key, min, max, step) => {
      const pair = document.createElement("span");
      pair.className = "dgxe-pair";
      const input = Object.assign(document.createElement("input"),
                                  {type: "range", min, max, step});
      input.style.accentColor = spec.colour;
      const out = document.createElement("output");
      const text = document.createElement("label");
      text.textContent = label;
      const set = (v) => {
        input.value = v;
        out.textContent = (+v).toFixed(step < 1 ? 1 : 0) + " \u00B0C";
      };
      input.addEventListener("input", () => {
        state[i][key] = +input.value;
        set(input.value);
        box.value = state.map((s) => ({...s}));
        box.dispatchEvent(new Event("input", {bubbles: true}));
      });
      pair.append(text, input, out);
      row.append(pair);
      return () => set(spec[key]);
    };
    row.prepend(name);
    const resets = [add("optimum", "optimum", 8, 24, 0.5),
                    add("upper limit", "limit", 28, 41, 1)];
    resets.forEach((r) => r());
    box.append(row);
    return resets;
  });
  const reset = document.createElement("button");
  reset.type = "button";
  reset.className = "dgxe-reset";
  reset.textContent = "Reset";
  reset.addEventListener("click", () => {
    state.forEach((s, i) => Object.assign(s, {optimum: specs[i].optimum,
                                              limit: specs[i].limit}));
    rows.forEach((r) => r.forEach((f) => f()));
    box.value = state.map((s) => ({...s}));
    box.dispatchEvent(new Event("input", {bubbles: true}));
  });
  box.append(reset);
  box.value = state.map((s) => ({...s}));
  return box;
}

// Remove the copy buttons of hidden code cells (outside the language tabs).
window.addEventListener("load", () => {
  for (const box of document.querySelectorAll(".code-copy-outer-scaffold")) {
    const pre = box.querySelector("pre");
    if (!box.closest(".tab-pane") && !(pre && pre.offsetHeight)) box.remove();
  }
});
