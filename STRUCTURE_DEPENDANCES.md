# Dependency Structure - Complete Organization

## Dependency Directories

### 📁 `include/utils/` - File Reading Utilities

The original dependencies have been moved to the appropriate directory:

```
include/utils/
├── readinteger.h          # Positive integer reading
├── readreal.h             # Positive real number reading
├── readnegativereal.h     # Real number reading (with sign and scientific notation)
├── FileIOUtils.h          # Convenience header (includes everything)
└── README.md              # Utilities documentation
```

## Complete Modular Organization

```
.
├── include/
│   ├── utils/                    # ⭐ NEW: Utility dependencies
│   │   ├── readinteger.h
│   │   ├── readreal.h
│   │   ├── readnegativereal.h
│   │   ├── FileIOUtils.h
│   │   └── README.md
│   │
│   ├── Constants.h
│   ├── ErrorCodes.h
│   ├── Exceptions.h
│   ├── Interfaces.h
│   ├── ChromosomeDivider.h
│   ├── GenomeDataManager.h
│   ├── RelativeIdentificationEngine.h
│   ├── GenomeFileLoader.h
│   ├── PhasingAlgorithmEngine.h
│   ├── OutputFileWriter.h
│   ├── ConfigurationManager.h
│   ├── HaplotypePhasingProgram.h
│   └── PhasingProgram.h
│
└── src/
    ├── GenomeDataManager.cpp
    ├── RelativeIdentificationEngine.cpp
    ├── GenomeFileLoader.cpp
    ├── ConfigurationManager.cpp
    ├── OutputFileWriter.cpp
    ├── PhasingAlgorithmEngine.cpp
    └── HaplotypePhasingProgram.cpp
```

## Dependency Usage

### Option 1: Convenience Header
```cpp
#include "include/utils/FileIOUtils.h"
// Automatically includes readinteger.h, readreal.h, readnegativereal.h
```

### Option 2: Individual Inclusion
```cpp
#include "include/utils/readinteger.h"
#include "include/utils/readreal.h"
#include "include/utils/readnegativereal.h"
```

### Option 3: Via PhasingProgram.h
```cpp
#include "include/PhasingProgram.h"
// Automatically includes FileIOUtils.h
```

## Reference Updates

All files have been updated to use the new paths:

- ✅ `GenomeFileLoader.h` → includes `utils/readinteger.h`
- ✅ `PhasingAlgorithmEngine.cpp` → includes `utils/readinteger.h`
- ✅ `PhasingProgram.h` → includes `utils/FileIOUtils.h`
- ✅ `Makefile_modular` → paths updated

## Advantages of this Organization

1. **Clear separation**: Utilities separated from main modules
2. **Reusability**: Utilities usable in other projects
3. **Documentation**: README.md dedicated to utilities
4. **Convenience header**: `FileIOUtils.h` for easy inclusion
5. **Logical organization**: All I/O utilities in the same place

## Migration

The old files at the root have been saved with the `.old` extension:
- `readinteger.h.old`
- `readreal.h.old`
- `readnegativereal.h.old`

They can be deleted after verifying that everything works correctly.
