#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

// Константы
#define NUM_CORNERS 4
#define EPSILON 1e-4

double euclidean_distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
// Функции для частных производных
double dphi_dx(double x, double y) {
    return 0.0;
}

double dphi_dy(double x, double y) {
    double numerator = 1 - 2*y;

    double denumerator = pow(2, 0.5)* pow(y-pow(y,2), 0.5);

    return numerator/denumerator; 
}

double dpsi_dx(double x, double y) {
    return -1/(2*pow(1-x, 0.5));
}

double dpsi_dy(double x, double y) {
    return 0;
}

// Частные производные для Df
double df_dx(double x, double y) {
    double numerator = 1.0; 
    double denumerator = 2*pow(1-x, 0.5); 
    return numerator/denumerator;
}

double df_dy(double x, double y) {
    return -1;
}

double dg_dx(double x, double y) {
    return 2*x;
}

double dg_dy(double x, double y) {
    return 4*y-2;
}

// Обратный Якобиан для Df
void inverse_jacobian(double x, double y, double J_inv[2][2]) {
    double df_dx_val = df_dx(x, y);
    double df_dy_val = df_dy(x, y);
    double dg_dx_val = dg_dx(x, y);
    double dg_dy_val = dg_dy(x, y);

    double det = df_dx_val * dg_dy_val - df_dy_val * dg_dx_val;
    if (fabs(det) < 1e-10) {
        printf("Ошибка: Якобиан необратим.\n");
        J_inv[0][0] = J_inv[1][1] = J_inv[0][1] = J_inv[1][0] = NAN;
        return;
    }

    J_inv[0][0] = dg_dy_val / det;
    J_inv[0][1] = -df_dy_val / det;
    J_inv[1][0] = -dg_dx_val / det;
    J_inv[1][1] = df_dx_val / det;
}

// Вычисление нормы матрицы (максимум суммы модулей элементов по строкам)
double matrix_norm(double J[2][2]) {
    double norm1 = fabs(J[0][0]) + fabs(J[0][1]);
    double norm2 = fabs(J[1][0]) + fabs(J[1][1]);
    return fmax(norm1, norm2);
}

// Функции для итерационного процесса
double phi_x(double x, double y) {
    return pow(2,0.5)*pow(y-pow(y,2),0.5);
}

double phi_y(double x, double y) {
    return pow(1-x, 0.5);
}

// Сами функции f и g
double f(double x, double y) {
    return pow(x+1, 0.5) - y;
}

double g(double x, double y) {
    return pow(x,2)+2*pow(y,2)-2*y;
}

// Вычисление q и mu для сегмента
void compute_q_mu_segments(const double segments[][4], int num_segments, double* q, double* mu) {
    for (int seg = 0; seg < num_segments; seg++) {
        double x_start = segments[seg][0];
        double x_end = segments[seg][1];
        double y_start = segments[seg][2];
        double y_end = segments[seg][3];

        double X[NUM_CORNERS] = { x_start, x_end, x_start, x_end };
        double Y[NUM_CORNERS] = { y_start, y_start, y_end, y_end };

        double norm_Dphi[NUM_CORNERS] = { 0 };
        double norm_Df[NUM_CORNERS] = { 0 };
        double norm_Df_inv[NUM_CORNERS] = { 0 };

        for (int i = 0; i < NUM_CORNERS; i++) {
            double x = X[i];
            double y = Y[i];

            double Jphi[2][2] = { {dphi_dx(x, y), dphi_dy(x, y)},
                                 {dpsi_dx(x, y), dpsi_dy(x, y)} };

            double Jf[2][2] = { {df_dx(x, y), df_dy(x, y)},
                               {dg_dx(x, y), dg_dy(x, y)} };

            double Jf_inv[2][2];
            inverse_jacobian(x, y, Jf_inv);

            norm_Dphi[i] = matrix_norm(Jphi);
            norm_Df[i] = matrix_norm(Jf);
            if (!isnan(Jf_inv[0][0])) {
                norm_Df_inv[i] = matrix_norm(Jf_inv);
            }
            else {
                norm_Df_inv[i] = NAN;
            }
        }

        double max_norm_Dphi = 0.0, max_norm_Df = 0.0, max_norm_Df_inv = 0.0;
        for (int i = 0; i < NUM_CORNERS; i++) {
            if (norm_Dphi[i] > max_norm_Dphi) max_norm_Dphi = norm_Dphi[i];
            if (norm_Df[i] > max_norm_Df) max_norm_Df = norm_Df[i];
            if (!isnan(norm_Df_inv[i]) && norm_Df_inv[i] > max_norm_Df_inv) max_norm_Df_inv = norm_Df_inv[i];
        }

        q[seg] = max_norm_Dphi;
        if (max_norm_Df_inv > 0) {
            mu[seg] = max_norm_Df * max_norm_Df_inv;
        }
        else {
            mu[seg] = NAN;
        }
    }
}

void seidel_method(double initial_points[1][2], double epsilon) {
    double x_prev = initial_points[0][0], y_prev = initial_points[0][1];
    double x_new, y_new;
    int iter = 0;

    while (true) {
        x_new = phi_x(x_prev, y_prev);
        y_new = phi_y(x_new, y_prev);

        //printf("x = %.4f, y = %.4f", x_new, y_new);

        if (euclidean_distance(x_prev, y_prev, x_new, y_new) < epsilon) {
            break;
        }

        x_prev = x_new;
        y_prev = y_new;
        iter++;
    }

    printf("Seidel method:\n");
    printf("root: x = %.4f, y = %.4f\n", x_new, y_new);
    printf("Iteration: %d\n", iter);
    printf("---------------------------\n");
}

// Итерационный метод
void iterative_method(double initial_points[1][2], double q, double epsilon) {
    double x_prev = initial_points[0][0], y_prev = initial_points[0][1];
    double x_new, y_new;
    int iter = 0;

    while (true) {
        x_new = phi_x(x_prev, y_prev);
        y_new = phi_y(x_prev, y_prev);

        //printf("x = %.4f, y = %.4f", x_new, y_new);

        if (euclidean_distance(x_prev, y_prev, x_new, y_new) < ((1 - q) / q) * epsilon) {
            break;
        }

        x_prev = x_new;
        y_prev = y_new;
        iter++;
    }

    printf("Iterative method:\n");
    printf("root: x = %.4f, y = %.4f\n", x_new, y_new);
    printf("Iteration: %d\n", iter);
    printf("---------------------------\n");
}

// Метод Ньютона
void newton_method(double initial_points[1][2], double mu, double epsilon) {
    double x_prev = initial_points[0][0], y_prev = initial_points[0][1];
    double x_new, y_new;
    int iter = 0;

    while (true) {
        double Jf_inv[2][2];
        inverse_jacobian(x_prev, y_prev, Jf_inv);

        double delta_x = -Jf_inv[0][0] * f(x_prev, y_prev) - Jf_inv[0][1] * g(x_prev, y_prev);
        double delta_y = -Jf_inv[1][0] * f(x_prev, y_prev) - Jf_inv[1][1] * g(x_prev, y_prev);

        x_new = x_prev + delta_x;
        y_new = y_prev + delta_y;

        if (euclidean_distance(x_prev, y_prev, x_new, y_new) < epsilon / mu) {
            break;
        }

        x_prev = x_new;
        y_prev = y_new;
        iter++;
    }

    printf("Newton method:\n");
    printf("root: x = %.4f, y = %.4f\n", x_new, y_new);
    printf("Iteration: %d\n", iter);
    printf("---------------------------\n");
}

int main() {
    double initial_points[1][2] = { { 0.66, 0.5 } };   // Начальные значения для сегмента 
    double segments[1][4] = { {0.5, 0.75, 0.25, 0.75} };  // Для одного сегмента
    int num_segments = 1;  // Количество сегментов

    // Динамическое выделение памяти под массивы q и mu
    double* q = (double*)malloc(num_segments * sizeof(double));
    double* mu = (double*)malloc(num_segments * sizeof(double));

    if (q == NULL || mu == NULL) {
        printf("Ошибка выделения памяти!\n");
        return 1;  // Возвращаем ошибку, если память не была выделена
    }

    // Вычисляем q и mu для всех сегментов
    compute_q_mu_segments(segments, num_segments, q, mu);

    // Применяем методы для каждого сегмента
    for (int i = 0; i < num_segments; i++) {
        printf("Segment %d: q = %.4f, mu = %.4f\n", i + 1, q[i], mu[i]);
        printf("---------------------------\n");

        // Выполняем методы для каждого сегмента
        seidel_method(initial_points, EPSILON);
        iterative_method(initial_points, 0.5, EPSILON);
        newton_method(initial_points, mu[i], EPSILON);
    }

    // Освобождение выделенной памяти
    free(q);
    free(mu);

    return 0;
}