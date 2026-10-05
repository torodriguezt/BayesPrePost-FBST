# Datos de entrada

`vaping_mccauley2023.csv` contiene conteos agregados y porcentajes publicados de
McCauley, Baiocchi, Cruse y Halpern-Felsher (2023). No son registros
individuales de estudiantes. La columna `source` registra la procedencia y
`status` distingue los conteos publicados de los reconstruidos.

`sample` distingue `linked` (estudiantes vinculados) y `full` (todos los
observados en cada etapa); `item` identifica la pregunta; `n1`, `x1`, `n2`, `x2`
son los márgenes. Cuando están disponibles, `n00`, `n01`, `n10` y `n11` dan la
tabla de transición. Los campos vacíos no se sustituyen por ceros.

La lectura y las comprobaciones están en
[`../Data_preparation/vaping_data.R`](../Data_preparation/vaping_data.R).
