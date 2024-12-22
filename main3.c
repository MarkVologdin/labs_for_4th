#include <stdio.h>
#include <math.h>
#include <locale.h>
#include <string.h>

double calculate_machine_epsilon() {
    double epsilon = 1.0;
    while ((1.0 + epsilon / 2.0) > 1.0) {
        epsilon /= 2.0;
    }
    return epsilon;
}

double third_derivative(const char* func_name, double x) {
    if (strcmp(func_name, "sin(x)") == 0) {
        return fabs(sin(x));
    }
    else if (strcmp(func_name, "10^5*sin(x)") == 0) {
        return fabs(100000 * sin(x));
    }
    else if (strcmp(func_name, "tan(x)") == 0) {
        return fabs(2*sin(x)/pow(cos(x), 3));
    }
    else {
        return 0.0;
    }
}

double check_error(const char* func_name, double x, double h) {
    double approx_derivative;
    if (strcmp(func_name, "sin(x)") == 0) {
        approx_derivative = (1 * sin(x) - 2 * sin(x + h) + sin(x + 2 * h)) / ( h*h);
        return fabs((approx_derivative + sin(x)) / sin(x));
    }
    else if (strcmp(func_name, "10^5*sin(x)") == 0) {
        approx_derivative = (-100000 * sin(x) - 200000 * sin(x + h) + 100000 * sin(x + 2 * h)) / (h * h);
        return fabs((approx_derivative + 100000 * sin(x)) / (100000 * sin(x)));
    }
    else if (strcmp(func_name, "tan(x)") == 0) {
        approx_derivative = (1 * tan(x) - 2 * tan(x + h) + tan(x + 2 * h)) / (h * h);
        return fabs((approx_derivative + 2*sin(x)/pow(cos(x),3) / (2*sin(x)/pow(cos(x),3))));
    }
    else {
        return 0.0;
    }
}

int main() {
    //setlocale(LC_ALL, "English");

    FILE* file = fopen("output.txt", "w"); // Open file for writing
    if (file == NULL) {
        printf("File opening error!\n");
        return 1; // Exit the program if the file could not be opened
    }

    double epsilon_machine = calculate_machine_epsilon();
    fprintf(file, "Machine epsilon (double): %.16e\n", epsilon_machine); // Write to file

    double x0 = 1.59; // Evaluation point
    double initial_h_min = pow(epsilon_machine, 1.0 / 3.0);
    double h_min; // Current step size
    double error;
    int k = 6; // Number of step size divisions
    double M3, M0;
    const char* functions[] = { "sin(x)", "10^5*sin(x)", "tan(x)" };
    int num_functions = sizeof(functions) / sizeof(functions[0]);

    for (int i = 0; i < num_functions; i++) {
        fprintf(file, "Function: %s\n", functions[i]);
        M3 = third_derivative(functions[i], x0);
        h_min = initial_h_min;
        for (int j = 0; j <= k; j++) {
            error = check_error(functions[i], x0, h_min);

            fprintf(file, "Step size h: %.10e\n", h_min);
            fprintf(file, "Approximation error for h: %.10e\n\n", error); // Write errors to file

            h_min /= 10.0;
        }
        fprintf(file, "====================================\n");
    }

    fclose(file); // Close file
    return 0;
}

