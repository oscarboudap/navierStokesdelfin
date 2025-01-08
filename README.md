# Simulación del Flujo de Fluidos alrededor de un Delfín

Este proyecto utiliza una aproximación numérica basada en las ecuaciones de Navier-Stokes para simular el flujo de fluidos alrededor de un delfín. Se genera una animación en formato GIF que muestra el comportamiento del flujo mientras interactúa con un obstáculo circular que representa al delfín.

## Características

- **Aproximación Numérica**: Utiliza diferencias finitas para resolver las ecuaciones de Navier-Stokes.
- **Obstáculo**: Representación del delfín como una máscara circular en el dominio.
- **Flujo Inicial**: Flujo uniforme inicial en la dirección \(x\).
- **Animación**: Se genera un GIF mostrando la evolución del flujo alrededor del delfín.

## Requisitos

Para ejecutar este proyecto necesitas tener instalado:

- **Python 3.8+**
- Bibliotecas:
  - `numpy`
  - `matplotlib`
  - `imagemagick` (para guardar el GIF)

### Instalación de dependencias
Ejecuta el siguiente comando para instalar las dependencias:
```bash
pip install numpy matplotlib
