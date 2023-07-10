# SLAP
Simple Linear Algebra Package (SLAP)

## Features
- Written in `C` for maximum compatibility among various systems
- No external libraries (all-in-one header file)
- Small enough to fit inside **MS-DOS** (and eventually embedded systems)
- Tailored to be used in _Finite Element_ software


# ToDo
- simple script that buld an ASCII file assessing the internal header-files structure (do inside the `header_merger`?)
- rendere coerenti i nomi delle funzioni (ad esempio mettendo `matd_` prima di ogni operazione sulle matrici double e poi il nome della funzione. Come ad esempio `matd_new`, `matd_free`, `matd_equal`, etc.)
- dividere la cartella `tests` da quella `examples`
- immagine di presentazione 1280×640px fatta carina con `SIL` che possa rendere accattivante cliccare su SLAP
- esempi nel README.md per come fare delle operazioni base con la libreria (come ad esempio risolvere un piccolo sistema lineare o qualche operazione base sui vettori)
