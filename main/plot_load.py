import matplotlib.pyplot as plt
import numpy as np
import datetime
import os



# Función para leer datos de un archivo de texto
def read_data(file_path):
    x = []
    y = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        if len(lines[0].split()) == 1:
            lines = lines[1:]
        for line in lines:
            data = line.split()
            x.append(float(data[0]))
            y.append(float(data[1]))
    return x, y

# Función para filtrar datos dentro de un rango específico
def filter_data(x, y, lower_bound=-160, upper_bound=-140):
    filtered_x = []
    filtered_y = []
    for i in range(len(y)):
        if lower_bound <= y[i] <= upper_bound:
            filtered_x.append(x[i])
            filtered_y.append(y[i])
    return filtered_x, filtered_y

# Función para graficar datos y guardar la imagen
def plot_data(x, y=None, title='Plot'):
    if y is None:
        y = x
        x = range(len(y))  # Generar valores x como un rango si no se proporcionan


    # Generar una marca de tiempo y crear un nombre de archivo
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    #filename = f"{title.replace(' ', '_')}_{timestamp}.png"
    filename = f"{title.replace(' ', '_')}.png"


    plt.figure()
    plt.plot(x, y, marker='o')
    if title != 'Plot':
        plt.title(title)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    #plt.title('Plot from text file')
    plt.grid(True)
    #plt.show()
    plt.savefig(filename)  # Guardar la figura como un archivo de imagen
    plt.close()
    os.startfile(filename)  # Abrir el archivo de imagen guardado con el visor de imágenes predeterminado



# Función para modificar los datos x en el archivo
# Borrar datos de primera columna y reemplazarlos por números pares
def modify_x_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    new_lines = []
    for i, line in enumerate(lines):
        data = line.split()
        data[0] = str(i * 2)
        new_lines.append(' '.join(data) + '\n')
    
    with open(file_path, 'w') as file:
        file.writelines(new_lines)

# Función para calcular el ruido usando un promedio móvil
def calculate_noise(x, y, window_size=100):
    y = np.array(y)
    moving_avg = np.convolve(y, np.ones(window_size)/window_size, mode='valid')
    noise = y[window_size-1:] - moving_avg
    x_adjusted = np.array(x[window_size-1:])/2
    return x_adjusted, noise


# Función para generar ruido adaptativo
def adaptative_noise(data, input_value, window_size=400, alpha=0.8):
    # Verificar que el tamaño de los datos sea mayor o igual al tamaño de la ventana
    if len(data) < window_size:
        raise ValueError("Data size is smaller than the window size.")
    
    # Seleccionar aleatoriamente una ventana de datos
    random_start = np.random.randint(0, len(data) - window_size)
    window = data[random_start:random_start + window_size]

    # Calcular el patrón de ruido restando la media de la ventana
    window_mean = np.mean(window)
    noise_pattern = window - window_mean

    # Normalizar el patrón de ruido
    max_abs_noise = np.max(np.abs(noise_pattern))
    normalized_noise = noise_pattern / (max_abs_noise + 1e-8)

    # Calcular la relación relativa entre el valor de entrada y la media de la ventana
    relative_ratio = np.abs(input_value / window_mean) if window_mean != 0 else 1

    # Calcular las distancias absolutas y los factores exponenciales
    distances = np.abs(noise_pattern)
    max_distance = np.max(distances)
    exponential_factors = np.exp(-distances / (max_distance + 1e-8))

    # Calcular los factores de escala combinando la relación relativa y los factores exponenciales
    scale_factors = alpha * relative_ratio + (1 - alpha) * exponential_factors * relative_ratio

    # Escalar y desnormalizar el ruido
    scaled_noise = normalized_noise * scale_factors * max_abs_noise
    denormalized_noise = scaled_noise + input_value

    # Graficar el ruido desnormalizado
    #plot_data(denormalized_noise, title='Denormalized Noise')

    return denormalized_noise


# Función para escribir un array en un archivo de texto con números pares en la primera columna
def write_array_to_file(array, file_path):
    with open(file_path, 'w') as file:
        for i, value in enumerate(array):
            file.write(f"{i * 2} {value}\n")


# Función principal
def main():
    #path = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\Code\main'  # Replace with your file path
    path = r'C:\Users\lldie\OneDrive - Universidad Técnica Federico Santa María\Universidad\Memoria\AGC\Escenarios\Salida ERNC'  # Replace with your file path

    name_file = r'\CHCs.txt'

    path_full = path + r'\Demanda.txt'
    path_short = path + r'\DatosDda\150.txt'
    #plot_data(*read_data(path_full), title='Full Data')
    #modify_x_data(path + r'\DatosDda\150.txt')

    plot_data(*read_data(path + name_file), title='Salida_CHCs')


    # Generar ruido adaptativo y guardarlo en un archivo
    value = 110
    y = adaptative_noise(read_data(path_full)[1], value)
    write_array_to_file(y, path + '\\DatosDda\\' + str(value) + '_noise.txt')
    plot_data(y, title='Adaptative Noise' + str(value))
    #plot_data(*read_data(path_short))
    

    #plot_data(*read_data(path_full))
    #x, y = filter_data(*read_data(path_full))
    #plot_data(x, y)

if __name__ == "__main__":
    main()
