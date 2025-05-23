from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import sys
import imageio  # для создания gif
import os
from ufl import And

# === Подготовка ===
mesh = Mesh("plate.xml")
boundaries = MeshFunction("size_t", mesh, "plate_facet_region.xml")
ds = Measure("ds", subdomain_data=boundaries)

deg = 2
V = FunctionSpace(mesh, "CG", deg)

# === Информации о сетке и числе искомых величин ===
n_c = mesh.num_cells()
n_v = mesh.num_vertices()
n_d = V.dim()
n = FacetNormal(mesh)  

print(f"Число ячеек сетки: {mesh.num_cells()}")
print(f"Число узлов сетки: {mesh.num_vertices()}")
print(f"Число искомых дискретных значений: {V.dim()}")

# === Пробная и тестовая функции ===
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)  

# === Граничные условия ===
u_l = Expression("x[1]", degree=deg)
u_r = Expression("x[1]", degree=deg)
u_t = Expression("x[1]", degree=deg)
u_b = Expression("0.0", degree=deg)

bcs1 = [#DirichletBC(V, u_l, boundaries, 2),
        #DirichletBC(V, u_r, boundaries, 3),
        DirichletBC(V, u_t, boundaries, 4),
        DirichletBC(V, u_b, boundaries, 5)]

# === Начальные условия ===
Gamma = -2.0
omega_0 = '-2.0'
omega_k = Function(V)
omega_k.interpolate(Expression(omega_0, degree=deg))
u_k = Function(V)

# === Подготовка визуализации ===
coordinates = mesh.coordinates()
x = coordinates[:, 0]
y = coordinates[:, 1]
triangles = mesh.cells()

image_files = []
output_dir = "frames"
os.makedirs(output_dir, exist_ok=True)

# === Итерационный процесс ===
max_iter = 10
tolerance = 1e-12

for k in range(max_iter):
    print(f"\n== Итерация {k} ==")
    a = dot(grad(u), grad(TestFunction(V))) * dx
    L = omega_k * TestFunction(V) * dx
    solve(a == L, u_k, bcs1)

    X, Y = SpatialCoordinate(mesh)
    indicator = conditional(
        lt(u_k, 0.0),
        conditional(And(gt(X, 1.65), And(gt(Y, 0.0), lt(Y, 0.5))),
                    1.0, 0.0),
        0.0
    )
    area_negative_u = assemble(indicator * dx)

    if near(area_negative_u, 0.0):
        print("Область ψ < 0 исчезла, остановка итераций")

        # === Визуализация и сохранение кадров ===
        u_values = u_k.compute_vertex_values(mesh)
        plt.figure(figsize=(15, 10))

        # Цветная подложка
        contourf = plt.tricontourf(x, y, triangles, u_values, levels=50, cmap='viridis')

        # Чёрные линии уровня
        contour = plt.tricontour(x, y, triangles, u_values,
                                levels=np.linspace(-1, 1, 40),
                                colors='black', linewidths=0.5)

        # Подписи значений на линиях
        plt.clabel(contour, fmt="%.2f", colors='black', fontsize=8)
        # Подписи осей
        plt.xlabel("L")
        plt.ylabel("H")

        # Цветовая шкала с подписью
        cbar = plt.colorbar(contourf)
        cbar.set_label("Функция тока ")
        plt.title(f"Iteration {k}")
        plt.colorbar(contourf)
        filename = f"{output_dir}/frame_{k:03d}.png"
        plt.savefig(filename)
        plt.close()
        image_files.append(filename)
        if near(area_negative_u, 0.0):
            print("Указанная область пуста.")
        else:
            print(f"Макс. завихренность: {np.max(omega_array):.6f}")
            print(f"Мин. завихренность: {np.min(omega_array):.6f}")
        break
        

    omega_expr = conditional(
        lt(u_k, 0.0),
        conditional(And(gt(X, 1.65), And(gt(Y, 0.0), lt(Y, 0.5))),
                    Constant(Gamma / area_negative_u),
                    Constant(0.0)),
        Constant(0.0)
    )

    omega_k = project(omega_expr, V)
    #omega_new = project(conditional(lt(u_k, 0.0), Constant(Gamma/area_negative_u), Constant(0.0)), V)
    omega_array = omega_k.vector().get_local()
    if near(area_negative_u, 0.0):
        print("Указанная область пуста.")
    else:
        print(f"Макс. завихренность: {np.max(omega_array):.6f}")
        print(f"Мин. завихренность: {np.min(omega_array):.6f}")

    # === Визуализация и сохранение кадров ===
    u_values = u_k.compute_vertex_values(mesh)
    plt.figure(figsize=(15, 10))

    # Цветная подложка
    contourf = plt.tricontourf(x, y, triangles, u_values, levels=50, cmap='viridis')

    # Чёрные линии уровня
    contour = plt.tricontour(x, y, triangles, u_values,
                            levels=np.linspace(-1, 1, 40),
                            colors='black', linewidths=0.5)

    # Подписи значений на линиях
    plt.clabel(contour, fmt="%.2f", colors='black', fontsize=8)

    # Подписи осей
    plt.xlabel("L")
    plt.ylabel("H")

    # Цветовая шкала с подписью
    cbar = plt.colorbar(contourf)
    cbar.set_label("Функция тока ")

    # Заголовок
    plt.title(f"Iteration {k}")

    # Сохранение
    filename = f"{output_dir}/frame_{k:03d}.png"
    plt.savefig(filename)
    plt.close()
    image_files.append(filename)

    # Проверка сходимости
    if k > 0:
        change = errornorm(u_k, u_prev, 'L2')
        print(f"Change: {change}")
        if change < tolerance:
            print("Достигнута сходимость")
            break
    u_prev = Function(V)
    u_prev.assign(u_k)


