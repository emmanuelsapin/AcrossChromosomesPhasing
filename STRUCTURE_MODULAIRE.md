# Modular Structure - Overview

## Summary

The code has been **split into a multitude of modular dependencies** for a professional and maintainable architecture.

## Dependency Architecture

### 📁 File Structure

```
.
├── include/                          # 13 header files
│   ├── Constants.h                  # System constants
│   ├── ErrorCodes.h                 # Error codes
│   ├── Exceptions.h                 # Custom exceptions
│   ├── Interfaces.h                 # Abstract interfaces (3 interfaces)
│   ├── ChromosomeDivider.h          # Data structure
│   ├── GenomeDataManager.h          # Genomic data management
│   ├── RelativeIdentificationEngine.h  # Relative identification
│   ├── GenomeFileLoader.h           # File loading
│   ├── PhasingAlgorithmEngine.h     # Phasing engine
│   ├── OutputFileWriter.h           # Output writing
│   ├── ConfigurationManager.h       # Configuration management
│   ├── HaplotypePhasingProgram.h    # Main program
│   └── PhasingProgram.h             # Main header (includes everything)
│
└── src/                             # 7 implementation files
    ├── GenomeDataManager.cpp
    ├── RelativeIdentificationEngine.cpp
    ├── GenomeFileLoader.cpp
    ├── ConfigurationManager.cpp
    ├── OutputFileWriter.cpp
    ├── PhasingAlgorithmEngine.cpp
    └── HaplotypePhasingProgram.cpp
```

## Dependency Hierarchy

### Level 1 - Foundations (No dependencies)
- ✅ `Constants.h`
- ✅ `ErrorCodes.h`

### Level 2 - Infrastructure
- ✅ `Exceptions.h` → depends on `ErrorCodes.h`
- ✅ `Interfaces.h` → pure interfaces (no dependencies)

### Level 3 - Data Structures
- ✅ `ChromosomeDivider.h` → simple structure (no dependencies)

### Level 4 - Data Management
- ✅ `GenomeDataManager.h` → depends on `Interfaces.h`, `Constants.h`
- ✅ `RelativeIdentificationEngine.h` → depends on `Interfaces.h`, `Constants.h`

### Level 5 - I/O
- ✅ `GenomeFileLoader.h` → no external dependencies
- ✅ `OutputFileWriter.h` → depends on `GenomeDataManager.h` (forward declaration)

### Level 6 - Algorithms
- ✅ `PhasingAlgorithmEngine.h` → depends on multiple modules
  - `Interfaces.h`
  - `ChromosomeDivider.h`
  - `Constants.h`
  - Forward declarations: `GenomeDataManager`, `RelativeIdentificationEngine`

### Level 7 - Configuration
- ✅ `ConfigurationManager.h` → no external dependencies

### Level 8 - Orchestration
- ✅ `HaplotypePhasingProgram.h` → depends on `ConfigurationManager.h`
  - Forward declarations for all other modules

### Level 9 - Entry Point
- ✅ `PhasingProgram.h` → convenience header including all modules

## Advantages of this Modular Architecture

### ✅ Modularity
- Each module can be compiled independently
- Unit tests per module
- Reuse in other projects

### ✅ Maintainability
- Isolated modifications to one module
- Incremental compilation
- Clear separation of responsibilities

### ✅ Extensibility
- Easy addition of new modules
- Module replacement without affecting others
- Support for multiple implementations (Strategy Pattern)

### ✅ Testability
- Unit tests per module
- Easy mocking via interfaces
- Modular integration tests

## Compilation

```bash
# Modular compilation
make -f Makefile_modular

# Clean
make -f Makefile_modular clean
```

## Usage

```cpp
// Using a specific module
#include "include/GenomeDataManager.h"
#include "include/RelativeIdentificationEngine.h"

using namespace PhasingEngine;

// Or complete usage
#include "include/PhasingProgram.h"
```

## Statistics

- **13 modular header files**
- **7 separate implementation files**
- **3 abstract interfaces** for decoupling
- **6 dependency levels** hierarchical
- **0 circular dependencies** (clean architecture)

## Next Steps

To add a new module:

1. Create `include/NewModule.h`
2. Create `src/NewModule.cpp`
3. Add to `Makefile_modular`
4. Update `PhasingProgram.h` if necessary
