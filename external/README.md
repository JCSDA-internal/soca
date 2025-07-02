# External Dependencies

This directory contains external submodules used by SOCA. Currently, SOCA is the only repository in JEDI that uses these dependencies (MOM6 and Icepack), which is why they are maintained as submodules here rather than as separate forks.

## Future Improvements

In the future, CMake logic should be added to handle these dependencies more efficiently when SOCA is built as part of a larger project:

- When built as part of a larger project, SOCA should link to existing MOM6 and ICEPACK installations rather than rebuilding them
- This would prevent duplicate builds and improve maintenance
