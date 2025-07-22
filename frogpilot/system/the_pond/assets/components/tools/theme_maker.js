import { html, reactive } from "https://esm.sh/@arrow-js/core";

export function ThemeMaker() {
  const state = reactive({
    themeName: "",
    turnSignalLength: 100,
    turnSignalType: "Traditional",
    colors: {
      LaneLines: { red: 23, green: 134, blue: 68, alpha: 242 },
      LeadMarker: { red: 23, green: 134, blue: 68, alpha: 242 },
      Path: { red: 23, green: 134, blue: 68, alpha: 242 },
      PathEdge: { red: 18, green: 107, blue: 54, alpha: 242 },
      Sidebar1: { red: 23, green: 134, blue: 68, alpha: 242 },
      Sidebar2: { red: 23, green: 134, blue: 68, alpha: 242 },
      Sidebar3: { red: 23, green: 134, blue: 68, alpha: 242 },
    },
    images: {
      homeButton: null,
      settingsButton: null,
      steeringWheel: null,
      turnSignal: null,
    },
    imageFileNames: {
      homeButton: "",
      settingsButton: "",
      steeringWheel: "",
      turnSignal: "",
    },
    sounds: {
      startup: null,
      prompt_repeat: null,
      engage: null,
      disengage: null,
    },
    soundFileNames: {
      startup: "",
      prompt_repeat: "",
      engage: "",
      disengage: "",
    },
  });

  function handleThemeNameInput(event) {
    state.themeName = event.target.value;
  }

  function handleImageUpload(event, key) {
    const file = event.target.files[0];
    if (file && !file.type.startsWith('image/')) {
      showSnackbar("Invalid file type! Please upload an image.", "error");
      event.target.value = "";
      return;
    }
    state.images[key] = file;
    state.imageFileNames[key] = file ? file.name : "";
  }

  function handleSoundUpload(event, key) {
    const file = event.target.files[0];
    if (file && !file.type.startsWith('audio/')) {
      showSnackbar("Invalid file type! Please upload an audio file.", "error");
      event.target.value = "";
      return;
    }
    state.sounds[key] = file;
    state.soundFileNames[key] = file ? file.name : "";
  }

  function handleColorChange(event, key) {
    const hex = event.target.value;
    const red = parseInt(hex.slice(1, 3), 16);
    const green = parseInt(hex.slice(3, 5), 16);
    const blue = parseInt(hex.slice(5, 7), 16);
    state.colors[key] = { red, green, blue, alpha: 255 };
  }

  function handleTurnSignalLengthInput(event) {
    state.turnSignalLength = event.target.value;
  }

  function validateTurnSignalLength(event) {
    let value = parseInt(event.target.value, 10);
    if (isNaN(value) || value < 25) {
      value = 25;
    } else if (value > 1000) {
      value = 1000;
    }
    state.turnSignalLength = value;
    event.target.value = value;
  }

  function toggleTurnSignalType(type) {
    state.turnSignalType = type;
  }

  async function saveTheme() {
    if (!state.themeName.trim()) {
      showSnackbar("Please enter a theme name...", "error");
      return;
    }

    const formData = new FormData();
    formData.append("themeName", state.themeName);
    formData.append("colors", JSON.stringify(state.colors));
    formData.append("turnSignalLength", state.turnSignalLength);
    formData.append("turnSignalType", state.turnSignalType);

    for (const key in state.images) {
      if (state.images[key]) {
        formData.append(key, state.images[key]);
      }
    }

    for (const key in state.sounds) {
      if (state.sounds[key]) {
        formData.append(key, state.sounds[key]);
      }
    }

    try {
      const response = await fetch("/api/themes", {
        method: "POST",
        body: formData,
      });

      const result = await response.json();

      if (response.ok) {
        showSnackbar(result.message || "Theme saved successfully!", "success");
      } else {
        showSnackbar(result.message || "Failed to save theme...", "error");
      }
    } catch (error) {
      showSnackbar("An error occurred while saving the theme...", "error");
    }
  }

  async function downloadTheme() {
    if (!state.themeName.trim()) {
      showSnackbar("Please enter a theme name to download.", "error");
      return;
    }

    const formData = new FormData();
    formData.append("themeName", state.themeName);
    formData.append("colors", JSON.stringify(state.colors));
    formData.append("turnSignalLength", state.turnSignalLength);
    formData.append("turnSignalType", state.turnSignalType);

    for (const key in state.images) {
      if (state.images[key]) {
        formData.append(key, state.images[key]);
      }
    }

    for (const key in state.sounds) {
      if (state.sounds[key]) {
        formData.append(key, state.sounds[key]);
      }
    }

    try {
      const response = await fetch("/api/themes/download", {
        method: "POST",
        body: formData,
      });

      if (response.ok) {
        const blob = await response.blob();
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `${state.themeName.replace(/ /g, "_")}.zip`;
        document.body.appendChild(a);
        a.click();
        a.remove();
        window.URL.revokeObjectURL(url);
      } else {
        const result = await response.json();
        showSnackbar(result.message || "Failed to download theme...", "error");
      }
    } catch (error) {
      showSnackbar("An error occurred while downloading the theme...", "error");
    }
  }

  async function deleteTheme() {
    if (!state.themeName.trim()) {
      showSnackbar("Please enter a theme name to delete.", "error");
      return;
    }

    if (confirm(`Are you sure you want to delete the theme "${state.themeName}"?`)) {
      try {
        const response = await fetch(`/api/themes/delete/${state.themeName}`, {
          method: "DELETE",
        });
        const result = await response.json();
        if (response.ok) {
          showSnackbar(result.message || "Theme deleted successfully!", "success");
          state.themeName = "";
        } else {
          showSnackbar(result.message || "Failed to delete theme.", "error");
        }
      } catch (error) {
        showSnackbar("An error occurred while deleting the theme.", "error");
      }
    }
  }

  return html`
    <div class="theme-maker-container">
      <div class="theme-maker-main-widget">
        <div class="theme-maker-main-title">Theme Maker</div>
        <div class="theme-name-section">
          <label for="themeName" class="theme-name-label">Theme Name</label>
          <input
            type="text"
            id="themeName"
            placeholder="Enter theme name..."
            autocomplete="off"
            value="${() => state.themeName}"
            @input="${handleThemeNameInput}"
          />
        </div>
        <div class="theme-maker-sub-widgets">
          <section class="theme-maker-widget">
            <div class="theme-maker-title">Colors</div>
            <div class="theme-maker-form">
              <div class="color-section">
                ${Object.keys(state.colors).sort().map(key => {
                  const labelText = {
                    LaneLines: "Lane Lines",
                    LeadMarker: "Lead Marker",
                    Path: "Path",
                    PathEdge: "Path Edge",
                    Sidebar1: "Sidebar Top",
                    Sidebar2: "Sidebar Middle",
                    Sidebar3: "Sidebar Bottom",
                  }[key];

                  return html`
                    <label class="color-label">
                      ${labelText}
                      <input
                        type="color"
                        value="${() => {
                          const color = state.colors[key];
                          return '#'
                            + color.red.toString(16).padStart(2, '0')
                            + color.green.toString(16).padStart(2, '0')
                            + color.blue.toString(16).padStart(2, '0');
                        }}"
                        @input="${e => handleColorChange(e, key)}"
                      />
                    </label>
                  `;
                })}
              </div>
            </div>
          </section>

          <section class="theme-maker-widget">
            <div class="theme-maker-title">Icons</div>
            <div class="theme-maker-form">
              <div class="upload-section">
                ${Object.keys(state.images).map(key => html`
                  <label class="file-upload-label">
                    <span class="file-upload-text">${{
                      homeButton: "Home Button",
                      settingsButton: "Settings Button",
                      steeringWheel: "Steering Wheel",
                      turnSignal: "Turn Signal",
                    }[key]}</span>
                    <span class="file-name-display">${() => state.imageFileNames[key]}</span>
                    <span class="file-upload-button">Choose File</span>
                    <input type="file" class="file-upload-input" accept="image/*"
                      @change="${e => handleImageUpload(e, key)}" />
                  </label>
                `)}
                <div class="turn-signal-length-section">
                  <label for="turnSignalLength" class="theme-name-label turn-signal-label">Turn Signal Length (25-1000ms)</label>
                  <input type="text" pattern="\\d*" id="turnSignalLength"
                    value="${() => state.turnSignalLength}"
                    @input="${handleTurnSignalLengthInput}"
                    @blur="${validateTurnSignalLength}"
                    class="turn-signal-input"
                  >
                </div>
                <div class="turn-signal-style-section">
                  <label class="theme-name-label turn-signal-label">Turn Signal Style</label>
                  <div class="signal-type-toggle">
                    <button
                      class="${() => `toggle-button ${state.turnSignalType === 'Static' ? 'active' : ''}`}"
                      @click="${() => toggleTurnSignalType('Static')}"
                    >Static</button>
                    <button
                      class="${() => `toggle-button ${state.turnSignalType === 'Traditional' ? 'active' : ''}`}"
                      @click="${() => toggleTurnSignalType('Traditional')}"
                    >Traditional</button>
                  </div>
                </div>
              </div>
            </div>
          </section>

          <section class="theme-maker-widget">
            <div class="theme-maker-title">Sounds</div>
            <div class="theme-maker-form">
              <div class="upload-section">
                ${[
                  { key: 'disengage', label: 'Disengage Sound' },
                  { key: 'engage', label: 'Engage Sound' },
                  { key: 'prompt_repeat', label: 'Prompt Sound' },
                  { key: 'startup', label: 'Startup Sound' }
                ].map(({ key, label }) => html`
                  <label class="file-upload-label">
                    <span class="file-upload-text">${label}</span>
                    <span class="file-name-display">${() => state.soundFileNames[key]}</span>
                    <span class="file-upload-button">Choose File</span>
                    <input type="file" class="file-upload-input" accept="audio/*"
                      @change="${e => handleSoundUpload(e, key)}" />
                  </label>
                `)}
              </div>
            </div>
          </section>
        </div>
        <div class="save-button-wrapper">
          <button class="delete-button" @click="${deleteTheme}">Delete Theme</button>
          <button class="download-button" @click="${downloadTheme}">Download Theme</button>
          <button class="save-button" @click="${saveTheme}">Save Theme</button>
        </div>
      </div>
    </div>
  `;
}
