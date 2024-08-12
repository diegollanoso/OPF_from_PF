# Code

### TBD

- Que tomar como perdidas? Solo líneas o tambien los trafos?? Se consideran ambos! No se pueden añadir perdidas en simulaciones 
- Rampas no sirven para desconexión de carga!


- Añadir perturbaciones de demanda.
- Eliminar carpetas y archivos innecesarios.

### Considerar

- Existen diferencias entre el SF obtenido con los flujos y algebraicamente. Posiblemente el error se de en el obtenido con el flujo, notar que aparece sensibilidad de potencia hacia el Gen 7_1

- Diferencias en las cargas, al final del CPF y se obtiene la demanda de las cargas durante la simulación, no se colocan condiciones. Mientras, que en el diccionario original se colocan 2 resitriccion, Out of Service / Circuit Breaker.

### Progreso


- 11/8: Se soluciono el problema de rampas
- 11/8: Se añadieron Excels con los resultados.

- 06/8: Añadir variación de demanda y perdidas por la variación de tensión. 
    - No olvidar la restricción de la potencia ENS <= DemandaPFC


- 05/8: Se añadio la variación de demanda por la diferencia de tensión en las cargas.

- 04/8: Se termino de añadir el efecto de disminución de potencia en unidades participantes del AGC de la formulación.

- 30/7: Se comenzo a añadir el efecto de disminución de potencia en unidades participantes del AGC de la formulación.