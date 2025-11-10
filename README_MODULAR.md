# Modular Architecture - Haplotype Phasing Engine

## Dependency Structure

The code has been split into a multitude of independent modules for better maintainability and reusability.

### Directory Structure

```
.
├── include/              # Header files (interfaces)
│   ├── Constants.h      # System constants
│   ├── ErrorCodes.h     # Error codes
│   ├── Exceptions.h     # Exception classes
│   ├── Interfaces.h     # Abstract interfaces
│   ├── ChromosomeDivider.h
│   ├── GenomeDataManager.h
│   ├── RelativeIdentificationEngine.h
│   ├── GenomeFileLoader.h
│   ├── PhasingAlgorithmEngine.h
│   ├── OutputFileWriter.h
│   ├── ConfigurationManager.h
│   ├── HaplotypePhasingProgram.h
│   └── PhasingProgram.h  # Main header (includes everything)
│
├── src/                 # Implementations
│   ├── GenomeDataManager.cpp
│   ├── RelativeIdentificationEngine.cpp
│   ├── GenomeFileLoader.cpp
│   ├── ConfigurationManager.cpp
│   └── OutputFileWriter.cpp
│
├── core/                 # Core modules (reserved for extensions)
├── data/                 # Data management (reserved)
├── io/                   # I/O operations (reserved)
└── config/              # Configuration (reserved)
```

## Module Dependencies

### Level 1 - Foundations
- **Constants.h**: No dependencies
- **ErrorCodes.h**: No dependencies
- **Exceptions.h**: Depends on `ErrorCodes.h`

### Level 2 - Interfaces
- **Interfaces.h**: No dependencies (pure interfaces)
- **ChromosomeDivider.h**: No dependencies (simple structure)

### Level 3 - Data Management
- **GenomeDataManager.h**: 
  - Depends on `Interfaces.h`
  - Depends on `Constants.h`
  - Depends on `Exceptions.h` (via implementation)

- **RelativeIdentificationEngine.h**:
  - Depends on `Interfaces.h`
  - Depends on `Constants.h`

### Level 4 - I/O and Algorithms
- **GenomeFileLoader.h**: No external dependencies
- **OutputFileWriter.h**: Depends on `GenomeDataManager.h` (forward declaration)
- **PhasingAlgorithmEngine.h**:
  - Depends on `Interfaces.h`
  - Depends on `ChromosomeDivider.h`
  - Depends on `Constants.h`
  - Depends on `GenomeDataManager.h` (forward declaration)
  - Depends on `RelativeIdentificationEngine.h` (forward declaration)

### Level 5 - Configuration and Orchestration
- **ConfigurationManager.h**: No external dependencies
- **HaplotypePhasingProgram.h**:
  - Depends on `ConfigurationManager.h`
  - Depends on all other modules (forward declarations)

### Level 6 - Entry Point
- **PhasingProgram.h**: Includes all headers (convenience header)

## Advantages of this Architecture

1. **Modularity**: Each module can be compiled and tested independently
2. **Reusability**: Modules can be reused in other projects
3. **Maintainability**: Isolated modifications to one module
4. **Testability**: Unit tests per module
5. **Incremental Compilation**: Only modified modules are recompiled
6. **Separation of Responsibilities**: Each module has a clear role

## Compilation

```bash
make -f Makefile_modular
```

## Adding New Modules

To add a new module:

1. Create `include/NewModule.h` with declarations
2. Create `src/NewModule.cpp` with implementations
3. Add dependencies in `PhasingProgram.h` if necessary
4. Update `Makefile_modular` with new files

## Usage Example

```cpp
#include "include/GenomeDataManager.h"
#include "include/RelativeIdentificationEngine.h"

using namespace PhasingEngine;

// Independent module usage
GenomeDataManager genomeManager;
RelativeIdentificationEngine relativeEngine;
```
