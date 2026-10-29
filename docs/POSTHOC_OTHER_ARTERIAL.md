# Análisis post hoc de otras arterias en Chile

Extensión solicitada el 07/10/2026, incorporada al paquete el 18/10/2026.
No corresponde a un análisis preespecificado. Los scripts anteriores y las
definiciones de los contrastes principales se mantienen.

## Ejecución

Desde la raíz del paquete:

```bash
python run_pipeline.py
```

Para ejecutar solamente la extensión:

```bash
python scripts/08_other_arterial_posthoc.py
```

Para validar todas las salidas ya generadas:

```bash
python scripts/06_validate.py
```

## Definiciones

- Componentes: I70/I73, I71, I72, I74 e I77–I79.
- Comparadores de A22: I71, otras arterias excluyendo I71, todo I70–I79 y coronario I20–I25.
- Venoso: I26 + I80–I82, sin I83.
- Período: 2015–2025; siete grupos de edad de 20–29 a ≥80.
- A21: razones de tasas agrupadas mujeres/hombres con denominadores INE; sin ajuste anual ni IC de RR.
- A22: OR caso–caso por edad e IC log-Wald; resumen MH por edad e IC de Robins–Breslow–Greenland.

Las tablas CSV conservan precisión numérica completa; redondear solo al
presentarlas. A21 tiene formato ancho; A22 tiene una fila por comparador
y grupo de edad o resumen MH. Los CSV de auditoría incluyen los conteos
y denominadores usados.

## Verificación

Se comprueba la partición de los 87.899 egresos I70–I79, los 109.820 egresos
venosos, y la concordancia por año, sexo y edad con el insumo principal.
Se cotejan los OR y los IC con los contrastes coronario y de otras arterias
ya presentes en el pipeline, incluidos los resúmenes MH. Los resultados
deben reproducir los valores de A21/A22 del manuscrito a su precisión
publicada. La referencia se conserva separada del cálculo con el nombre
y SHA-256 del manuscrito fuente. El informe completo está en
`results/validation.json`.

Esta comprobación no valida clínicamente los códigos ni identifica personas
únicas. Las tasas son tasas de egreso. Los intervalos del OR son nominales
y no corrigen reingresos correlacionados. Las causas del patrón por edad
no pueden inferirse de esta descomposición.

## Repositorio

El paquete está listo para incorporarse a un repositorio. Esta actualización
local no publica ni sube archivos a GitHub. Se conservan los manifiestos
de fuentes históricas; el manifiesto PACKAGE_SHA256.csv se actualiza para
identificar los bytes de la nueva entrega.
