import streamlit as st
from PIL import Image
import base64
import subprocess
import os
import math
import gmsh
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

def run_gmsh(file_path):
    try:
        env = os.environ.copy()
        env["LIBGL_ALWAYS_SOFTWARE"] = "1"  # Используем программный рендеринг
        subprocess.run(["gmsh", file_path], check=True, env=env)
        st.success("Gmsh успешно запущен в программном режиме!")
    except FileNotFoundError:
        st.error("Gmsh не найден. Убедитесь, что он установлен и доступен в PATH.")
    except subprocess.CalledProcessError:
        st.error("Ошибка при запуске Gmsh.")

def show_code(code, language="python"):
    st.code(code, language)

menu = st.sidebar.radio(
    '***',
    (
        "Расчетная область",
        "Пример",
    )
)

if menu == "Расчетная область":

    r"""
    ##### Геометрическая модель
    **Расчетная область**

    Симметричное обтекание

        """

    c1, c2 = st.columns([3,1])
    c1.image("11.png")

    r"""
        
        + $h - $ полувысота пластины
        + $(l_1, l_2) - $ положение пластины   
        + $(l_2 - l_1) - $ ширина пластины   
        + $\Omega - $ расчетная область  
        + $H, L - $ полуширина канала, расчетная длина канала
        + $\Gamma_1 - $ линии симметрии и граница пластины
        + $\Gamma_2 - $ твердые стенки канала   
        + $\Gamma_3, \Gamma_4 - $ граница входа и выхода жидкости

    """
if menu ==   "Пример":

    r"""
    ##### Файл геометрии
    """
    
    with st.expander("Параметры геометрии"):
        code_1 = """
        // Параметры прямоугольника
        l1 = 1.45;        // Длина до пластины
        l2 = 0.1;        // Ширина пластины
        l3 = 4.0 - l2 - l1; // Длина после пластины
        H  = 1.0;       // Высота канала
        h_plate = 0.5;   // Высота вертикальной пластины
        d = 0.01;        // Характерный размер сетки
        """
        show_code(code_1, "python")
    with st.expander("Общая длина"):
        code_2 = """
        L = l1 + l2 + l3;
        """
        show_code(code_2, "python")
    with st.expander("Построение точек канала"):
        code_3 = """
        // Размеры сетки
        Point(1) = {0, 0, 0, d};            // Нижний левый угол
        Point(2) = {L, 0, 0, d};            // Нижний правый угол
        Point(3) = {L, H, 0, d};            // Верхний правый угол
        Point(4) = {0, H, 0, d};            // Верхний левый угол
        """
        show_code(code_3, "python")
    with st.expander("Построение точек пластины"):
        code_4 = """
        Point(5) = {l1, 0, 0, d};
        Point(6) = {l1 + l2, 0, 0, d};
        Point(7) = {l1 + l2, h_plate, 0, d};
        Point(8) = {l1, h_plate, 0, d};
        """
        show_code(code_4, "python")
    with st.expander("Построение линий"):
        code_5 = """
        //Линии (нижняя грань с пластиной)
        Line(1) = {1, 5};        // нижняя левая до пластины
        Line(2) = {5, 8};        // левая сторона пластины
        Line(3) = {8, 7};        // верхняя сторона пластины
        Line(4) = {7, 6};        // правая сторона пластины
        Line(5) = {6, 2};        // нижняя правая после пластины
        Line(6) = {2, 3};        // правая стенка
        Line(7) = {3, 4};        // верхняя стенка
        Line(8) = {4, 1};        // левая стенка
        """
        show_code(code_5, "python")
    with st.expander("Построение контура и поверхности"):
        code_6 = """
        Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
        Plane Surface(1) = {1};
        """
        show_code(code_6, "python")
    with st.expander("Определение физических ГУ"):
        code_7 = """
         // Физические группы для FEniCS
        Physical Surface("fluid") = {1};

        Physical Curve("inlet")  = {8};
        Physical Curve("outlet") = {6};            
        Physical Curve("top") = {7};
        Physical Curve("bottom")  = {1, 2, 3, 4, 5};   
        """
        show_code(code_7, "python")
    with st.expander("Построение сетки"):
        code_8 = """
        // Настройка сетки
        Mesh.CharacteristicLengthMin = 0.1;
        Mesh.CharacteristicLengthMax = 0.25;
        Mesh 2;
        """
        show_code(code_8, "python")

    result = ''
    for i in range(1, 8):
        result += globals()[f'code_{i}']

    def save_example_file():
        example_file_path = './plate.geo'
        with open(example_file_path, 'w') as f:
                f.write(result)
        return example_file_path

    # Кнопка для загрузки и запуска примера
    if st.button("Построение расчетной области 🔧"):
            example_file_path = save_example_file()
            run_gmsh(example_file_path)
        
