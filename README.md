# Code

### TBD

- Realizar presentación de paper

- Adaptar la función objetivo y restricciones para añadir distintos escenarios de demanda.

- Añadir generadores estaticos.

- Necesito data de generación variable.

- Que tomar como perdidas? Solo líneas o tambien los trafos?? Se consideran ambos! No se pueden añadir perdidas en simulaciones 

- Adecuar restricciones de rampa cuando existe desconexión de carga.

- ¿las perdidas PL_PRE son obtenidas antes de realizar la contingencia?

- Idea: Solo extraer la información de los generadores que participan en el AGC.

- Añadir perturbaciones de demanda.
- Eliminar carpetas y archivos innecesarios.

### Considerar

- Se asume que todos las contigencias afectan al sistema, es decir, antes de la salida de cargas o generadores se encuentran consumiendo o entregando potencia.

- Añadir escenarios de demanda distintos al modelo (valle, promedio, peak)
    - Posibles soluciones: 
        1.- Usar listas para las cargas y generadores.
        2.- Exportar datos de potencia de cargas y generadores e importarlos como tuplas.

- Existen diferencias entre el SF obtenido con los flujos y algebraicamente. Posiblemente el error se de en el obtenido con el flujo, notar que aparece sensibilidad de potencia hacia el Gen 7_1

- Diferencias en las cargas, al final del CPF y se obtiene la demanda de las cargas durante la simulación, no se colocan condiciones. Mientras, que en el diccionario original se colocan 2 resitriccion, Out of Service / Circuit Breaker.

### Progreso

- 13/8: Se solucionan los SF, problemas en el control de frecuencia - potencia del generador. Se añaden los escenarios valle y mean.

- 12/8: Se obtiene SF de forma algebraica, existen diferencias con el obtenido por PowerFactory

- 11/8: Se soluciono el problema de rampas
- 11/8: Se añadieron Excels con los resultados.

- 06/8: Añadir variación de demanda y perdidas por la variación de tensión. 
    - No olvidar la restricción de la potencia ENS <= DemandaPFC


- 05/8: Se añadio la variación de demanda por la diferencia de tensión en las cargas.

- 04/8: Se termino de añadir el efecto de disminución de potencia en unidades participantes del AGC de la formulación.

- 30/7: Se comenzo a añadir el efecto de disminución de potencia en unidades participantes del AGC de la formulación.