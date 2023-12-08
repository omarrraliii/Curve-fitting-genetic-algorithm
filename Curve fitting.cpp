#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
using namespace std;
#define population_size 10
#define num_of_generations 100
#define PC 0.7
#define PM 0.1
#define k 2
// Define the fitness function to calculate the distance between points and curve
double MSE(const vector<double> chromosome, const vector<pair<double, double>> points, int num_points, int degree)
{
    double total_fitness = 0;
    for (int counter_of_points = 0; counter_of_points < num_points; counter_of_points++)
    {
        double fitness = 0;
        for (int i = 0; i < degree + 1; i++)
        {
            fitness = fitness + (chromosome[i] * pow(points[counter_of_points].first, i));
        }
        fitness = fitness - points[counter_of_points].second;
        fitness = fitness * fitness;
        total_fitness = total_fitness + fitness;
    }
    total_fitness = total_fitness / num_points;
    return total_fitness;
}

// Genetic Algorithm functions
vector<double> initialize_individual(int degree)
{
    vector<double> chromosome(degree + 1);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i < degree + 1; i++)
    {
        chromosome[i] = dist(gen);
    }
    return chromosome;
}
vector<vector<double>> initialize_population(int degree)
{
    vector<vector<double>> population(population_size);
    for (int i = 0; i < population_size; i++)
    {
        population[i] = initialize_individual(degree);
    }
    return population;
}

vector<vector<double>> crossover(const vector<double> parent1, const vector<double> parent2)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(1, static_cast<int>(parent1.size()) - 1);
    int crossover_point1 = dis(gen);
    int crossover_point2 = dis(gen);
    vector<double> child1(parent1);
    vector<double> child2(parent2);
    if (crossover_point1>crossover_point2)
        swap(crossover_point1,crossover_point2);
    for (int i = crossover_point1; i < crossover_point2; i++)
    {
        child1[i] = parent2[i];
    }
    for (int i = crossover_point1; i < crossover_point2; i++)
    {
        child2[i] = parent1[i];
    }
    return {child1, child2};
}

void mutate(vector<double> &individual, int gen)
{
    random_device rd;
    mt19937 regen(rd());
    uniform_real_distribution<double> dis(0.0, 0.999999999999);
    for (double &gene : individual)
    {
        if (dis(regen) <= PM)
        {
            double y = 0;
            double delta_L = gene + 10;
            double delta_U = 10 - gene;
            double random_number = dis(regen);
            if (random_number < 0.5)
                y = delta_L;
            else
                y = delta_U;
            double delta_T_Y = y * (1 - pow(dis(regen), (pow(1 - (gen / num_of_generations), 1))));
            if (random_number < 0.5)
                gene = gene - delta_T_Y;
            else
                gene = delta_T_Y - gene;
        }
    }
}
vector<vector<double>> elitisReplacement (vector<vector<double>> current_generation, vector<vector<double>> new_generation,
                                          vector<pair<double, int>> fitness_values_current_generation,vector<pair<double, int>> fitness_values_new_generation){
        sort (fitness_values_current_generation.rbegin(),fitness_values_current_generation.rend());
        sort (fitness_values_new_generation.rbegin(),fitness_values_new_generation.rend());
        vector<vector<double>> population(population_size);
        for (int i=0;i<k;i++){
            population[i]=current_generation[fitness_values_current_generation[i].second];
        }
        for (int i=k;i<population_size;i++){
            population[i]=new_generation[fitness_values_new_generation[i-k].second];
        }
        return population;
    }

vector<double> genetic_algorithm(const vector<pair<double, double>> &data_points, int degree, int num_points, double mutation_rate = 0.1)
{
    vector<vector<double>> population(population_size);
    population = initialize_population(degree);

    for (int gen = 0; gen < num_of_generations; ++gen)
    {
        vector<pair<double, int>> fitness_values_current_generation(population_size);
        vector<pair<double, int>> fitness_values_new_generation(population_size);
        for (int i = 0; i < population_size; ++i)
        {
            fitness_values_current_generation[i] = {1.0/(MSE(population[i], data_points, num_points, degree)), i};
        }
        random_device rd;                                                                       //tournament selection
        mt19937 regen(rd());
        uniform_int_distribution<int> dist(0, population_size-1);
        vector<pair<double, int>> mating_pool(population_size);
        for (int i = 0; i < population_size; i++)
        {
            int parent1_index = dist(regen);
            int parent2_index = dist(regen);
            if (fitness_values_current_generation[parent1_index].first >= fitness_values_current_generation[parent2_index].first)
            {
                mating_pool[i] = {fitness_values_current_generation[parent1_index].first, fitness_values_current_generation[parent1_index].second};
            }
            else
            {
                mating_pool[i] = {fitness_values_current_generation[parent2_index].first, fitness_values_current_generation[parent2_index].second};
            }
        }
        sort(mating_pool.rbegin(), mating_pool.rend());
        vector<vector<double>> offspring(population_size);
        for (int i = 0; i <= population_size - 2; i += 2)
        {
            vector<vector<double>> parents(2);
            parents[0] = population[mating_pool[i].second];
            parents[1] = population[mating_pool[i + 1].second]; 

            random_device rd;
            mt19937 regen(rd());
            uniform_real_distribution<double> dist(0.0, 0.99999999999999);
            vector<vector<double>> children(2);
            if (dist(regen) < PC)
            {
                children = crossover(parents[0], parents[1]);
                mutate(children[0], gen);
                mutate(children[1], gen);
            }
            else{
                children[0]=parents[0];
                children[1]=parents[1];
            }
            offspring[i] = children[0];
            offspring[i + 1] = children[1];
            fitness_values_new_generation[i] = {1.0/(MSE(population[i], data_points, num_points, degree)), i};
            fitness_values_new_generation[i+1] = {1.0/(MSE(population[i], data_points, num_points, degree)), i+1};
        }

        population = elitisReplacement(population,offspring,fitness_values_current_generation,fitness_values_new_generation);
    }
    

    // Find the best individual (fitting coefficients)
    vector<pair<double,int>>fitness_values_current_generation(population_size);
    for (int i = 0; i < population_size; ++i)
        {
            fitness_values_current_generation[i] = {1.0/(MSE(population[i], data_points, num_points, degree)), i};
        }
    sort (fitness_values_current_generation.rbegin(),fitness_values_current_generation.rend());
    vector<double> best_individual = population[fitness_values_current_generation[0].second];
    return best_individual;
}

int main()
{
    int num_datasets;
    ifstream inputFile("curve_fitting_input.txt");
    inputFile >> num_datasets;
    ofstream writer("curve_fitting_output.txt");
    for (int dataset = 0; dataset < num_datasets; ++dataset)
    {
        int num_points, degree;
        inputFile >> num_points >> degree;

        vector<pair<double, double>> data(num_points);
        for (int i = 0; i < num_points; ++i)
        {
            double x, y;
            inputFile >> x >> y;
            data[i] = {x, y};
        }

        // Apply genetic algorithm to find best coefficients for curve fitting
        std::vector<double> best_coefficients = genetic_algorithm(data, degree, num_points);
        
        writer << "Best coefficients for dataset " << dataset + 1 << ": ";
        for (double coeff : best_coefficients)
        {
            writer << coeff << " ";
        }
        writer << endl;
        writer << "MSE is: "<<MSE(best_coefficients,data,num_points,degree)<<endl;
    }
    writer.close();
    inputFile.close();
    return 0;
}