# Debugging Guide with VS Code + GDB

Welcome to the Ne2x Debugging Guide! This document will help you set up and use Visual Studio Code (VS Code) with GDB to debug the Ne2x program effectively.

## Prerequisites
- **VS Code**: Make sure you have [Visual Studio Code](https://code.visualstudio.com/) installed.
- **C/C++ Extension**: Install the official [C/C++ extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) by Microsoft.
- **GDB**: Ensure you have GDB installed on your system. You can check by running `gdb --version` in your terminal.
- **Ne2x Source Code**: You should have the Ne2x source code available on your machine.

## Debug Configurations

The `.vscode/launch.json` file contains several debug configurations tailored for Ne2x. Here are the main ones:

### 1. **Debug Ne2x (with i:info + o:option)** ⭐ Main
   - Executes: `./Ne2x i:test_info.txt o:test_option.txt`
   - Use this to debug normal mode with info + option files

### 2. **Debug Ne2x (interactive mode)**
   - Executes: `./Ne2x` (no arguments)
   - Use this to test interactive mode

### 3. **Debug Ne2x (multi files)**
   - Executes: `./Ne2x m:multi_input.txt`
   - Use this to test multiple-files mode

### 4. **Debug Ne2x (stop at main)** 🔍
   - Like #1 but stops at the beginning of main
   - Useful for stepping through from the start

## How to Use the Debugger

### Quick Start

1. Open `Ne2x.c` in the editor
2. Press `F5` or click "Run and Debug" in the sidebar
3. Select a configuration from the dropdown menu
4. The program will start with GDB active!

### Essential Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `F5` | Start/Continue debugging |
| `F9` | Toggle breakpoint (on current line) |
| `F10` | Step Over (execute line, don't enter functions) |
| `F11` | Step Into (enter functions) |
| `Shift+F11` | Step Out (exit current function) |
| `Ctrl+Shift+F5` | Restart debugging |
| `Shift+F5` | Stop debugging |

### Setting Breakpoints

**Method 1 - Click:**
- Click to the left of the line number (a red dot appears)

**Method 2 - Keyboard:**
- Position cursor on line and press `F9`

**Method 3 - Conditional:**
- Right-click on breakpoint → "Edit Breakpoint"
- Add condition, e.g.: `argc > 2` or `p == 5`

### Panels During Debug

When you start debugging, VS Code automatically opens:

1. **VARIABLES** (left sidebar)
   - See all local variables
   - Expand pointers and structs
   - Hover over variables in code to see them

2. **WATCH** (left sidebar)
   - Add expressions to monitor
   - E.g.: `argv[1]`, `FileOne`, `hasOpt`
   - Click "+" to add

3. **CALL STACK** (left sidebar)
   - See the call chain
   - Click to jump to a function in the stack

4. **BREAKPOINTS** (left sidebar)
   - Manage all breakpoints
   - Enable/disable without deleting them

5. **DEBUG CONSOLE** (bottom)
   - Execute GDB commands directly
   - Evaluate expressions
   - E.g.: `-exec print argc` or simply `argc`

## Practical Example: Tracing Your Command Flow

### Scenario: `./Ne2x i:test_info.txt o:test_option.txt`

1. **Set strategic breakpoints:**
   ```
   Line ~229: int main(int argc, char *argv[])
   Line ~290: if (argv[1][0] == 'i') {
   Line ~310: RunOption(...)
   ```

2. **Start debug**: Press `F5`, select "Debug Ne2x (with i:info + o:option)"

3. **When it stops at first breakpoint:**
   - VARIABLES panel shows: `argc = 3`, `argv = [...]`
   - In WATCH add: `argv[1]`, `argv[2]`
   - You'll see: `argv[1] = "i:test_info.txt"`, `argv[2] = "o:test_option.txt"`

4. **Continue (F5 or F10):**
   - You'll reach `if (argv[1][0] == 'i')`
   - See in VARIABLES that `FileOne` and `FileTwo` get populated

5. **Enter RunOption (F11):**
   - When you reach the `RunOption(...)` call, press `F11`
   - You'll enter inside the function
   - See all passed parameters

## Advanced Tips

### 1. Conditional Breakpoints
```
Right-click on breakpoint → Edit Breakpoint
Condition: argc == 3 && argv[1][0] == 'i'
```

### 2. Data Breakpoints (watchpoints)
```
Right-click on variable in VARIABLES panel → "Break on Value Change"
The debugger stops when the variable changes!
```

### 3. Logpoints (printf without modifying code)
```
Right-click to the left of line → "Add Logpoint"
Message: argc = {argc}, argv[1] = {argv[1]}
```

### 4. Debug Console Commands
In the DEBUG CONSOLE you can:
```
argc                    // See value
argv[1]                 // See string
*FileOne                // Dereference pointer
p argc                  // GDB print
info locals             // All local variables
backtrace               // Complete stack
```

### 5. Memory View (if you have C/C++ extension)
```
During debug: Right-click on pointer → "View Memory"
See the hexadecimal in memory!
```

## Modifying Program Arguments

To test with different arguments:

1. Open `.vscode/launch.json`
2. Find the configuration you use
3. Modify the `"args"` array:
   ```json
   "args": [
       "i:other_file.txt",
       "o:other_options.txt"
   ]
   ```
4. Or create a new configuration (copy-paste)

## Compile Flags Explained

In `tasks.json` the "Build Ne2x Debug" task uses:

- `-g3`: Maximum debug symbols (includes macros)
- `-O0`: No optimization (code 1:1 with source)
- `-DDEBUG`: Define DEBUG (you can use `#ifdef DEBUG`)
- `-Wall -Wextra`: All warnings

## Troubleshooting

### "Unable to start debugging"
→ Compile first: `Ctrl+Shift+B` or Terminal → Run Build Task

### Breakpoint "gray" (not active)
→ Program not compiled with `-g`, recompile in Debug mode

### "No symbol table is loaded"
→ Executable doesn't have debug symbols, use "Build Ne2x Debug" task

### Variables show `<optimized out>`
→ Compile with `-O0` (already done in Build Ne2x Debug)

## Agent Workflow for C/C++ Debugging

This project includes AI customization files that complement GDB debugging.

- Instruction file: `.github/instructions/c-memory-safety.instructions.md`
- Prompt file: `.github/prompts/run-core-validation.prompt.md`
- Skill file: `.github/skills/investigate-regression-mismatch/SKILL.md`

### When to use each one

1. `c-memory-safety.instructions.md`
   - Use when editing memory-sensitive C/C++ code (`Ne2x.c`, related headers, wrapper boundary).
   - It enforces safe ownership and cleanup patterns and reminds the required validation flow.

2. `/run-core-validation`
   - Use after any non-trivial C/C++ change.
   - It runs the standard sequence from repo root:
     - `make`
     - `make test`
     - `make asan`
     - `./Ne2x_asan i:test_info.txt o:test_option.txt`
   - It returns a compact PASS/FAIL report and the first failing excerpt.

3. `/investigate-regression-mismatch`
   - Use when `make test` fails or C and C++ outputs diverge.
   - It helps classify differences (timestamp-only vs scientific mismatch), point to likely code areas, and propose the smallest safe fix.

### Recommended debug loop

1. Reproduce with GDB using one of the launch configurations.
2. Fix the issue with minimal localized edit (prefer C core parity-safe changes).
3. Run `/run-core-validation`.
4. If regression still fails, run `/investigate-regression-mismatch`.
5. Re-run the same debug scenario to confirm behavior and stability.

## Quick Start for Today

```bash
# 1. Compile in debug mode (automatic on first F5)
# 2. Open Ne2x.c
# 3. Set a breakpoint on this line (around 229):
#    int main(int argc, char *argv[])
# 4. Press F5
# 5. Select "Debug Ne2x (with i:info + o:option)"
# 6. Explore with F10 (step over) and F11 (step into)!
```

---

Happy debugging.
