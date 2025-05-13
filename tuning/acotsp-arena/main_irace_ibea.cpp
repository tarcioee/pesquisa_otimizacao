#ifndef MAIN_CPP
#define MAIN_CPP

#include "main.h"
#include "BoundedParetoSet.cpp"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> 
#include <chrono>
#include <filesystem>

using namespace std;
using namespace std::chrono;

// basic functions

Instance get_instance(string file_path) {
  Instance instance;

  ifstream file(file_path.c_str());
  if (!file.is_open()) {
    cerr << "Erro ao abrir o arquivo " << file_path << endl;
    return instance;
  }

  // Lendo dados do arquivo
  file >> instance.number_of_cities >> instance.number_of_passengers >> instance.car_capacity;

  // Redimensionando matrizes
  instance.cost_matrix.resize(instance.number_of_cities,  vector<int>(instance.number_of_cities));
  instance.time_matrix.resize(instance.number_of_cities,  vector<int>(instance.number_of_cities));

  // Lendo cost_matrix
  for (int i = 0; i < instance.number_of_cities; ++i) {
    for (int j = 0; j < instance.number_of_cities; ++j) {
      file >> instance.cost_matrix[i][j];
    }
  }

  // Lendo time_matrix
  for (int i = 0; i < instance.number_of_cities; ++i) {
    for (int j = 0; j < instance.number_of_cities; ++j) {
      file >> instance.time_matrix[i][j];
    }
  }

  // Lendo passageiros
  int max_cost, origin, destiny, max_time;
  for (int i = 0; i < instance.number_of_passengers; i++) {
    Passenger passenger;
    file >> max_cost >> origin >> destiny >> max_time;   
    passenger.max_cost = max_cost;
    passenger.origin = origin;
    passenger.destiny = destiny;
    passenger.max_time = max_time;
    instance.passengers.push_back(passenger);
  }

  // Lendo min_quota
  file >> instance.min_quota;

  // Lendo bonus_and_time
  int city, bonus, time;
  while (file >> city >> bonus >> time) {
    instance.bonus_and_time.push_back(make_pair(bonus, time));
  }

  file.close();

  return instance;
} // extrai a instância de determinado file path e retornar uma instância


void print_fronts(Population population) {
  for(int i=0; i<population.fronts.size(); i++){
    cout << "Fronte " << i << ":" << endl;
    for(int j=0; j<population.fronts[i].size(); j++){
      cout<<population.fronts[i][j]<<" ";
    }
    cout<<endl;
  }
}

int random_chance() {
  int random_chance = rand() % 100;
  return random_chance;
}
// returns true if margin>random percentage and else otherwise

int random_city(Instance instance) {
  return rand() % instance.number_of_cities;
}
// retorna uma cidade aleatória da instância

void swap_random_route_slots(Solution &solution) {
  int i = rand() % solution.route.size();
  int j = rand() % solution.route.size();
  swap(solution.route[i], solution.route[j]);
}
// troca dois slots aleatórios de uma rota de uma solution

void swap(std::pair<double, int> &a, std::pair<double, int> &b) {
  pair<double, int> temp = a;
  a = b;
  b = temp;
}

void swap(int &a, int &b) {
  int temp = a;
  a = b;
  b = temp;
}

void swap(double &a, double &b) {
  double temp = a;
  a = b;
  b = temp;
}

// Função para ajustar o heap
void heapify(vector<pair<double, int>> &vec, int n, int i) {
  int largest = i;       // Inicializa o maior como raiz
  int left = 2 * i + 1;  // Esquerda = 2*i + 1
  int right = 2 * i + 2; // Direita = 2*i + 2

  // Se o filho esquerdo for maior que a raiz
  if (left < n && vec[left].first > vec[largest].first) {
    largest = left;
  }

  // Se o filho direito for maior que o maior até agora
  if (right < n && vec[right].first > vec[largest].first) {
    largest = right;
  }

  // Se o maior não é a raiz
  if (largest != i) {
    swap(vec[i], vec[largest]);
    // Recursivamente faz o heapify na subárvore afetada
    heapify(vec, n, largest);
  }
}

void heapify(vector<int> &vec, int n, int i) {
  int largest = i;       // Inicializa o maior como raiz
  int left = 2 * i + 1;  // Esquerda = 2*i + 1
  int right = 2 * i + 2; // Direita = 2*i + 2

  // Se o filho esquerdo for maior que a raiz
  if (left < n && vec[left] > vec[largest]) {
    largest = left;
  }

  // Se o filho direito for maior que o maior até agora
  if (right < n && vec[right] > vec[largest]) {
    largest = right;
  }

  // Se o maior não é a raiz
  if (largest != i) {
    swap(vec[i], vec[largest]);
    // Recursivamente faz o heapify na subárvore afetada
    heapify(vec, n, largest);
  }
}

void heapify(vector<double> &vec, int n, int i) {
  int largest = i;       // Inicializa o maior como raiz
  int left = 2 * i + 1;  // Esquerda = 2*i + 1
  int right = 2 * i + 2; // Direita = 2*i + 2

  // Se o filho esquerdo for maior que a raiz
  if (left < n && vec[left] > vec[largest]) {
    largest = left;
  }

  // Se o filho direito for maior que o maior até agora
  if (right < n && vec[right] > vec[largest]) {
    largest = right;
  }

  // Se o maior não é a raiz
  if (largest != i) {
    swap(vec[i], vec[largest]);
    // Recursivamente faz o heapify na subárvore afetada
    heapify(vec, n, largest);
  }
}

// Função principal para fazer o Heap Sort
void heapSort(vector<pair<double, int>> &vec) {
  int n = vec.size();

  // Constrói o heap (reorganiza o vetor)
  for (int i = n / 2 - 1; i >= 0; i--) {
    heapify(vec, n, i);
  }

  // Extrai um elemento do heap de cada vez
  for (int i = n - 1; i >= 0; i--) {
    // Move a raiz atual para o final
    swap(vec[0], vec[i]);
    // Chama o heapify na heap reduzida
    heapify(vec, i, 0);
  }
}// ORDENA CRESCENTE pelo first

void heapSort(vector<int> &vec) {
  int n = vec.size();

  // Constrói o heap (reorganiza o vetor)
  for (int i = n / 2 - 1; i >= 0; i--) {
    heapify(vec, n, i);
  }

  // Extrai um elemento do heap de cada vez
  for (int i = n - 1; i >= 0; i--) {
    // Move a raiz atual para o final
    swap(vec[0], vec[i]);
    // Chama o heapify na heap reduzida
    heapify(vec, i, 0);
  }
}

void heapSort(vector<double> &vec) {
  int n = vec.size();

  // Constrói o heap (reorganiza o vetor)
  for (int i = n / 2 - 1; i >= 0; i--) {
    heapify(vec, n, i);
  }

  // Extrai um elemento do heap de cada vez
  for (int i = n - 1; i >= 0; i--) {
    // Move a raiz atual para o final
    swap(vec[0], vec[i]);
    // Chama o heapify na heap reduzida
    heapify(vec, i, 0);
  }
}

bool x_dominates_y(Solution solutionx, Solution solutiony) {
  if (solutionx.cost <= solutiony.cost && solutionx.time <= solutiony.time &&
      solutionx.total_bonus >= solutiony.total_bonus) {
    if (solutionx.cost < solutiony.cost || solutionx.time < solutiony.time ||
        solutionx.total_bonus > solutiony.total_bonus) {
      return true;
    }
  }
  return false;
}

Population get_non_dominated_population(Population received_population){
  Population population = received_population;
  for(int solution = 0; solution < population.population.size(); solution++){
    for(int solution_to_compare = solution+1; solution_to_compare  < population.population.size(); solution_to_compare++){
      if(x_dominates_y(population.population[solution], population.population[solution_to_compare])){
        population.population.erase(population.population.begin() + solution_to_compare );
        solution_to_compare--;
      }
      else if(x_dominates_y(population.population[solution_to_compare], population.population[solution])){
        population.population.erase(population.population.begin() + solution );
        solution--;
        break;
      }
      else if(population.population[solution].cost == population.population[solution_to_compare].cost and population.population[solution].total_bonus == population.population[solution_to_compare].total_bonus and population.population[solution].time == population.population[solution_to_compare].time){
        population.population.erase(population.population.begin() + solution_to_compare );
        solution_to_compare--;
      }
    }
  }
  return population;
}

Pareto_objectives get_pareto_objectives(Population received_population){
  Population population = get_non_dominated_population(received_population);
  Pareto_objectives matrix_of_objectives;
  for(int solution =0;solution<population.population.size();solution++){
   Objectives vector_of_objectives;
   vector_of_objectives.cost = population.population[solution].cost;
   vector_of_objectives.time = population.population[solution].time;
   vector_of_objectives.total_bonus = population.population[solution].total_bonus;

   matrix_of_objectives.pareto_set.push_back(vector_of_objectives);
  }
  return matrix_of_objectives;
}


//versão que receberá BoundedParetoSet
Pareto_objectives get_pareto_objectives(BoundedParetoSet *received_population){
  Pareto_objectives matrix_of_objectives;
  int pop_size = received_population->get_size();
  for(int solution =0;solution<pop_size;solution++){
    matrix_of_objectives.pareto_set.push_back(received_population->get_objectives(solution));
    matrix_of_objectives.solutions.push_back(received_population->get_solution(solution));
  }
  return matrix_of_objectives;
}

int get_bonus(Instance instance, Solution &solution) {
  int bonus = 0;
  for (int city = 0; city < solution.route.size(); city++) {
    if (solution.cities_colected[city]) { // add bonus from cities collected
      bonus += instance.bonus_and_time[solution.route[city]].first;
    }
  }
  return bonus;
}
// retorna o bonus total da solução no .total_bonus

int get_time(Instance instance, Solution &solution) {
  int time = 0;
  for (int city = 0; city < solution.route.size(); city++) {
    int next_city;
    if(city==solution.route.size()-1){
      next_city=0;
    }
    else{
      next_city=city+1;
    }
    time += instance.time_matrix[solution.route[city]][solution.route[next_city]];

    if (solution.cities_colected[city]){ // add time from bonus collecting
      time += instance.bonus_and_time[solution.route[city]].second;
    }
  } 
  return time;
}
// retorna o tempo total da solução no .time

double get_cost(Instance instance, Solution &solution) {
  double cost = 0;
  int people_in_car = 1;

  for (int city = 0; city < solution.route.size(); city++) {
    double temp_cost = 0;
    int next_city;
    if(city==solution.route.size()-1){
      next_city=0;
    }
    else{
      next_city=city+1;
    }
    temp_cost = instance.cost_matrix[solution.route[city]][solution.route[next_city]];
    if(solution.passengers_riding.size()>0){
      for (int passenger = 0; passenger < instance.number_of_passengers; passenger++) {
        if (solution.passengers_riding[passenger]==true) {
          if (instance.passengers[passenger].origin == solution.route[city]) {
            people_in_car++;
          }
          if (instance.passengers[passenger].destiny == solution.route[city]) {
            if(city !=0){
              people_in_car--;
            }
          }
        }
      }
    }
    temp_cost /= people_in_car;
    cost += temp_cost;
    /*cout<<"temp_cost: "<<temp_cost<<" people in car: "<<people_in_car<<endl;*/
  } // at every city checks how many are in the car to calculate the cost for
    // this part of the route
  return cost;
}
// retorna o custo total da solução no .cost

void update_objectives(Instance instance, Solution &solution) {
  solution.time = get_time(instance, solution);
  solution.cost = get_cost(instance, solution);
  solution.total_bonus = get_bonus(instance, solution);
}
// chama as 3 funções de consultar valor dos objetivos e atualiza nos
// respectivos lugares

void greedy_route(Instance instance, Solution &solution) {
  vector<int> visited_cities;
  solution.route[0] = random_city(instance);
  int actualI = 0, actualj = 0, marker = 0;
  double value_holder;
  while (actualI < solution.route.size()) {
    for (int j = 1; j < solution.route.size(); j++) {
      if (visited_cities[j] != j) {
        if (marker == 0) {
          value_holder = instance.cost_matrix[solution.route[actualI]][j];
          marker++;
          actualj = j;
        } else if (instance.cost_matrix[solution.route[actualI]][j] <
                   value_holder) {
          value_holder = instance.cost_matrix[solution.route[actualI]][j];
          actualj = j;
        }
      }
    }
    actualI++;
    solution.route[actualI] = actualj;
    visited_cities[actualj] = actualj;
    marker = 0;
  }
}
// Constrói uma rota de forma gulosa respeitando o tamanho limite da rota

vector<int> get_random_route(Instance instance) {
  int number_of_cities = 0;
  while (number_of_cities < 2) {
    number_of_cities = rand() % instance.number_of_cities;
  } // define um numero de cidades >= 2
  vector<int> route;
  vector<int> visited_cities(instance.number_of_cities, -1);
  route.push_back(0); //satisfazer origem sempre 0
  visited_cities[0] =0; 
  int i = 1;
  while (i < number_of_cities) {
    int city = random_city(instance);
    if (visited_cities[city] != city) {
      route.push_back(city);
      visited_cities[city] = city;
      i++;
    }
  }
  return route;
}
// number of cities will be at minimum 2, max all cities, equal chances to all
// possibilities

vector<bool> get_random_bonus(Instance instance, Solution solution) {
  vector<bool> cities_colected;
  cities_colected.push_back(false);//primeira cidade não é coletada
  for (int i = 0; i < solution.route.size()-1; i++) {
    if (rand() % 2 == 0) {
      cities_colected.push_back(true);
    } else {
      cities_colected.push_back(false);
    }
  }
  return cities_colected;
}
// 50% de chance de coletar um bônus. Trata cities_colected como tendo mesmo
// tamanho de route. Isso implica que a cidade em route[0] foi coletada se
// cities_colected[0] for True.

double getObj(Solution s, int k){
	if(k == 0){
		return s.cost;
	}
	else if(k==1){
		return s.time;
	}
	else if(k==2){
		return s.total_bonus;
	}
  else{
    return -1;
  }
}

void print_solution(Solution solution) {
  cout << "Rota: ";
  for(int i = 0; i < solution.route.size(); i++){
    cout << solution.route[i] << " ";
  }
  cout<<endl;
  cout << "Custo: " << solution.cost << endl;
  cout << "Bonus: " << solution.total_bonus << endl;
  cout << "Tempo: " << solution.time << "\n \n";

}

void print_solution(Instance instance, Solution &solution) {
  cout << "Rota: ";
  for(int i = 0; i < solution.route.size(); i++){
    cout << solution.route[i] << " ";
  }
  cout<<endl;

  cout<<"custos da rota:"<<endl;
  for(int city = 0; city < solution.route.size(); city++){
    int next_city;
    if(city == solution.route.size() - 1){
      next_city = 0;
    }
    else{
      next_city = city+1;
    }
    cout << instance.cost_matrix[solution.route[city]][solution.route[next_city]]<< " ";
  }
  cout<<endl;

  cout<<"Tempos da rota:"<<endl;
  for(int city = 0; city < solution.route.size(); city++){
    int next_city;
    if(city == solution.route.size() - 1){
      next_city = 0;
    }
    else{
      next_city = city+1;
    }
    cout << instance.time_matrix[solution.route[city]][solution.route[next_city]]<< " ";
  }
  cout<<endl;

  cout << "Cities_colected: "<<endl;
  bool any_city_colected = false;
  for(int i =0; i<solution.cities_colected.size();i++){
    if(solution.cities_colected[i]==true){
      cout<<"City: "<< solution.route[i] << " bonus " << instance.bonus_and_time[solution.route[i]].first << " tempo " << instance.bonus_and_time[solution.route[i]].second << endl;
      any_city_colected = true;
    }
  }
  if(!any_city_colected){
    cout<<"Nenhuma cidade coletada"<<endl;
  }
  cout<<endl;
  int passengers_on = 0;
  for (int i = 0; i < solution.passengers_riding.size(); i++) {
    if (solution.passengers_riding[i]) {
      passengers_on++;
      cout<< "passageiro "<<i<<" com origem "<<instance.passengers[i].origin<<", destino: "<<instance.passengers[i].destiny<< " custo max: "<<instance.passengers[i].max_cost<<" tempo max: "<<instance.passengers[i].max_time<<endl;
    }
  }
  cout << "Passageiros embarcados: " << passengers_on << endl;
  cout << "Custo: " << solution.cost << endl;
  cout << "Bonus: " << solution.total_bonus << endl;
  cout << "Tempo: " << solution.time << "\n \n";
}

// Função para calcular o custo do passageiro
double calculate_passenger_cost(int origem_index, int destiny_index, vector<int> passengers_in_car_by_city, Solution solution, Instance instance) {
    double cost = 0.0;
    if (destiny_index == 0) {
        cost += instance.cost_matrix[solution.route[solution.route.size() - 1]][solution.route[0]] / (1+passengers_in_car_by_city[solution.route.size() - 1]);
        destiny_index = solution.route.size() - 1;
    }
    for (int i = origem_index; i < destiny_index; i++) {
        int next_index = i + 1;
        cost += instance.cost_matrix[solution.route[i]][solution.route[next_index]] / (1+ passengers_in_car_by_city[i]);
    }
    // Caso o destino seja a cidade inicial, incluir o custo do último trecho
    return cost;
}

// Função para calcular o tempo do passageiro
double calculate_passenger_time(int origem_index, int destiny_index, Solution solution, Instance instance) {
    double time = 0;
    // Caso o destino seja a cidade inicial, incluir o tempo do último trecho
    if (destiny_index == 0) {
        time += instance.time_matrix[solution.route[solution.route.size() - 1]][solution.route[0]];
        destiny_index = solution.route.size() - 1;
    }
    for (int i = origem_index; i < destiny_index; i++) {
        int next_index = i + 1;
        time += instance.time_matrix[solution.route[i]][solution.route[next_index]];
        if (solution.cities_colected[i]) { //embarque e desembarque é após coleta de bonus
            time += instance.bonus_and_time[solution.route[i]].second;
        }
    }
    if (solution.cities_colected[destiny_index]) { //embarque e desembarque é após coleta de bonus
            time += instance.bonus_and_time[solution.route[destiny_index]].second;
    }
    return time;
}

void able_passengers(Instance instance, Solution &solution) {
    vector<bool> able_passengers(instance.number_of_passengers, false); 
    vector<int> passengers_in_car_by_city(solution.route.size(), 1); 
    vector<pair<double, int>> passengers_by_cost;

    for(int i = 0; i < instance.number_of_passengers; i++) {
        passengers_by_cost.push_back(make_pair(instance.passengers[i].max_cost, i));
    }

    heapSort(passengers_by_cost);  // Ordena passageiros por custo (do mais pobre ao mais rico)

    for(int passenger = passengers_by_cost.size() - 1; passenger >= 0; passenger--) {  // Vai do mais rico pro mais pobre
        int origem_index = -1, destiny_index = -1;
        double cost = 0.0;
        double time = 0;
        bool subiu = false, desceu = false, apto = true;

        for (int city = 0; city < solution.route.size(); city++) {
            if (solution.route[city] == instance.passengers[passengers_by_cost[passenger].second].origin) {
                // Passageiro pode subir
                if (passengers_in_car_by_city[city] < instance.car_capacity) {
                    origem_index = city;
                    subiu = true;
                } else {
                    apto = false;
                    break;
                }
            } else if(solution.route[city] == instance.passengers[passengers_by_cost[passenger].second].destiny) {
                // Passageiro desce
                if(city==0){ //caso especial em que destino é cidade inicial
                    desceu = true;
                    destiny_index = city;
                }
                else{
                  if (subiu) { 
                      destiny_index = city;
                      desceu = true;
                      break;
                  } else {
                      apto = false;
                      break;
                  }
                }
            } else if (subiu) {
                if (passengers_in_car_by_city[city] >= instance.car_capacity) {// Se o carro estiver cheio em qualquer cidade no caminho
                  apto = false;
                  break;
                }//checar se passageiro faz parte do car capacity. se nao fizer, mudar ali pra >
            }
        }

        // Trata o caso em que o destino do passageiro é a cidade inicial (rota que retorna ao ponto de origem)
        if (subiu && desceu && solution.route[0] == instance.passengers[passengers_by_cost[passenger].second].destiny) {
            if(passengers_in_car_by_city[solution.route.size()-1] >= instance.car_capacity){ //testa se o carro está cheio
              apto = false;
            }
        }

        // FASE DE TESTES OBJETIVOS: Se o passageiro subiu e desceu, calcula o custo e tempo e testa se ele poder arcar
        if (subiu && desceu && apto) {
            cost = calculate_passenger_cost(origem_index, destiny_index, passengers_in_car_by_city, solution, instance);
            time = calculate_passenger_time(origem_index, destiny_index, solution, instance);

            // Verifica se o passageiro pode arcar com o custo e tempo
            if (cost > instance.passengers[passengers_by_cost[passenger].second].max_cost || 
                time > instance.passengers[passengers_by_cost[passenger].second].max_time) {
                apto = false;
            }
        }

        // FASE PÓS-TESTES: Se o passageiro for apto, marca como "able" e marca a ocupação do carro nos trechos em que estava
        if (subiu && desceu && apto) {
            able_passengers[passengers_by_cost[passenger].second] = true;
            if (destiny_index == 0) {
                passengers_in_car_by_city[solution.route.size()-1]++;
                destiny_index = solution.route.size() - 1;
            }
            for (int city = origem_index; city < destiny_index; city++){
                passengers_in_car_by_city[city]++;
            }
        }
    }
    solution.passengers_riding = able_passengers;
}



void able_passengers(Instance instance, Population &population) {
  for(int solution = 0; solution<population.population.size();solution++){
    able_passengers(instance, population.population[solution]);
  }
} // aplica heuristica de carregamento de passageiros em todos os individuos da população


bool check_passengers_riding(Instance instance, Solution solution) {
    vector<int> passengers_in_car_by_city(solution.route.size(), 1); 
    vector<pair<double, int>> passengers_by_cost;
    for(int i =0; i <instance.number_of_passengers;i++){//adiciona todos os passageiros carregados da solução em pares de custo e indice
      if(solution.passengers_riding[i]){ 
        passengers_by_cost.push_back(make_pair(instance.passengers[i].max_cost,i));
      }
    }
    if(passengers_by_cost.size()==0){
      return true;
    }
    heapSort(passengers_by_cost);//aqui já ordenamos os passageiros por mais pobre 

    for(int passenger = passengers_by_cost.size() - 1; passenger >= 0; passenger--) { // Vai do mais rico pro mais pobre
        int origem_index = -1, destiny_index = -1;
        double cost = 0.0;
        int time = 0;
        bool subiu = false, desceu = false, apto = true;

        for (int city = 0; city < solution.route.size(); city++) {
            if (solution.route[city] == instance.passengers[passengers_by_cost[passenger].second].origin) {
                // Passageiro pode subir
                if (passengers_in_car_by_city[city] < instance.car_capacity) {
                    origem_index = city;
                    subiu = true;
                } else {
                    return false;
                }
            }else if(solution.route[city] == instance.passengers[passengers_by_cost[passenger].second].destiny) {
                // Passageiro desce
                if(city==0){ //caso especial em que destino é cidade inicial
                    desceu = true;
                    destiny_index = city;
                }
                else{
                  if (subiu) { 
                      destiny_index = city;
                      desceu = true;
                      break;
                  } else {
                      return false;
                  }
                }
            }else if (subiu) {
                if (passengers_in_car_by_city[city] >= instance.car_capacity) {// Se o carro estiver cheio
                  return false;
                  break;
                }
            }
        }

        // Trata o caso em que o destino do passageiro é a cidade inicial (rota que retorna ao ponto de origem)
        if (subiu && desceu && solution.route[0] == instance.passengers[passengers_by_cost[passenger].second].destiny) {
            if(passengers_in_car_by_city[solution.route.size()-1] >= instance.car_capacity){ //testa se o carro está cheio
              return false;
            }
        }

        // FASE DE TESTES OBJETIVOS: Se o passageiro subiu e desceu, calcula o custo e tempo e testa se ele poder arcar
        if (subiu && desceu && apto) {
            cost = calculate_passenger_cost(origem_index, destiny_index, passengers_in_car_by_city, solution, instance);
            time = calculate_passenger_time(origem_index, destiny_index, solution, instance);

            // Verifica se o passageiro pode arcar com o custo e tempo
            if (cost > instance.passengers[passengers_by_cost[passenger].second].max_cost || 
                time > instance.passengers[passengers_by_cost[passenger].second].max_time) {
                return false;
            }
        }

        // FASE PÓS-TESTES: Se o passageiro for apto, marca a ocupação do carro nos trechos em que estava
        if (subiu && desceu && apto) {
            if (destiny_index == 0) {
                passengers_in_car_by_city[solution.route.size()-1]++;
                destiny_index = solution.route.size() - 1;
            }
            for (int city = origem_index; city < destiny_index; city++){
                passengers_in_car_by_city[city]++;
            }
        }
    }

    return true;
}


bool solution_validity(Instance instance, Solution solution){ //checar as outras restrições
  bool validity = true;
  if(solution.route.size()<2){
    validity = false;
    //cout<<"rota menor que 2"<<endl;
  }
  if(solution.route[0]!=0){
      validity = false;
      //cout<<"Origem não é 0"<<endl;
  }
  if(solution.cost==0 or solution.time ==0){
    validity = false;
    //cout<<"custo ou tempo 0"<<endl;
  }
  if(solution.cost != get_cost(instance, solution)){
    validity = false;
    //cout<<"custo errado"<<endl;
  }
  if(solution.time != get_time(instance, solution)){
    validity = false;
    //cout<<"tempo errado"<<endl;
  }
  if(solution.total_bonus != get_bonus(instance, solution)){
    validity = false;
    //cout<<"bonus errado"<<endl;
  }
  if(check_passengers_riding(instance, solution)==false){
    validity = false;
    //cout<<"passageiros errado"<<endl;
  }
    //if(validity ==false){ print_solution(instance,solution);}

  return validity;
}

Solution get_random_solution(Instance instance) {
  Solution solution;
  solution.route = get_random_route(instance);
  solution.cities_colected = get_random_bonus(instance, solution);
  able_passengers(instance, solution);
  update_objectives(instance, solution);
  return solution;
}

Population get_random_population(Instance instance, int max_population) {
  Population population;
  int i =0;
  while(i < max_population) {
    Solution solution;
    solution = get_random_solution(instance);
    if(solution_validity(instance, solution)){
      population.population.push_back(solution);
      i++;
    }
  }
  return population;
}

Population get_random_population(Instance instance, BoundedParetoSet *EP ,int max_population, int &valuations) {
  Population population;
  int i =0;
  while(i < max_population) {
    Solution solution;
    solution = get_random_solution(instance);
    valuations++;
    if(solution_validity(instance, solution)){
      population.population.push_back(solution);
      if(EP->adicionarSol(&solution)){
          continue;
      }
      i++;
    }
    
  }
  return population;
}

void mutate_routes(Instance instance, Solution &Solution, int mode) {
  if (mode == 3) {
    if (Solution.route.size() <= 2) {
      mode = 2;
    } else if (Solution.route.size() == instance.number_of_cities) {
      mode = rand() % 2;
    } else {
      mode = rand() % 3;
    }
  } // esse if assegura que se mutate_routes for chamado, com certeza ocorre
    // alguma mutação. talvez o professor prefira a possibilidade da chamada se
    // perder.
  if (mode == 0) { // swap vertex (but not origin)
    int city1 = rand() % Solution.route.size(); 
    while(city1 == 0){
      city1 = rand() % Solution.route.size();
    } 
    int city2 = rand() % Solution.route.size();
    while(city2==0){
      city2 = rand() % Solution.route.size();
    }
    swap(Solution.route[city1], Solution.route[city2]);
    swap(Solution.cities_colected[city1], Solution.cities_colected[city2]);
  } else if (mode == 1) { // remove random vertex, but not origin
      int city_to_diminish = 0;
      while(city_to_diminish==0){
        city_to_diminish = rand() % Solution.route.size();
      }
      Solution.route.erase(Solution.route.begin() + city_to_diminish);
      Solution.cities_colected.erase(Solution.cities_colected.begin() + city_to_diminish);
  }else if (mode == 2) { // add random vertex that is not already in the route in any place except origin
  
      vector<int> cities_not_in_route;
      bool marker = false;
      for (int j = 0; j < instance.number_of_cities; j++) {
        marker = false;
        for (int k = 0; k < Solution.route.size(); k++) {
          if (Solution.route[k] == j) {
            marker = true;
            break;
          }
        }
        if (marker == false) {
          cities_not_in_route.push_back(j);
        }
      }
      if (cities_not_in_route.size() == 0) {
        //cout << "deu errado pq tinha todos na rota: " << cities_not_in_route.size() << endl;
        //cout << "os da instancia: " << instance.number_of_cities << endl;
      }
      int city_to_add = rand() % cities_not_in_route.size();
      int city_to_add_index = rand() % Solution.route.size();
      while(city_to_add_index==0){
        city_to_add_index = rand() % Solution.route.size();
      }
      Solution.route.insert(Solution.route.begin() + city_to_add_index, cities_not_in_route[city_to_add]);
      Solution.cities_colected.insert(Solution.cities_colected.begin() + city_to_add_index, false);
  }
}
// mutates a solution in four different settings according to "mode": 0-
// swapping, 1 - removing, 2 - adding a new vertex and 3 - random

void mutate_routes(Instance instance, Population &population) {
  for (int i = 0; i < population.population.size(); i++) {
    mutate_routes(instance, population.population[i], 3);
  }
}
// mutates a population in three different ways: swapping, removing and adding a
// new vertex

void invert_bonuses(Instance instance, Solution &solution) {
  for (int i = 0; i < solution.cities_colected.size(); i++) {
    if (solution.cities_colected[i] == true) {
      solution.cities_colected[i] = false;
    } else {
      solution.cities_colected[i] = true;
    }
  }
}

void mutate_bonuses(Instance instance, Solution &solution) {
  int bonus_change = rand() % solution.cities_colected.size();
  while(bonus_change==0){
    bonus_change = rand() % solution.cities_colected.size();
  }
  if (solution.cities_colected[bonus_change] == false) {
    solution.cities_colected[bonus_change] = true;
  } else {
    solution.cities_colected[bonus_change] = false;
  }
}

void mutate_bonuses(Instance instance, Population &population) {
  for (int i = 0; i < population.population.size(); i++) {
    mutate_bonuses(instance, population.population[i]);
  }
}

Solution one_crossover(Instance instance, vector<Solution> &population, int father, int mother) {
  Solution baby;
  for (int i = 0; i < population[father].route.size() / 2; i++) {
    baby.route.push_back(population[father].route[i]);
    baby.cities_colected.push_back( population[father].cities_colected[i]);
  } // primeira metade da rota e coletas do pai

  for (int j = population[mother].route.size() / 2;
       j < population[mother].route.size(); j++) {
    bool gene_is_repeated = false;
    for (int k = 0; k < population[father].route.size() / 2; k++) {
      if (population[mother].route[j] ==  population[father].route[k]) {
        gene_is_repeated = true;
      }
    }
    if (gene_is_repeated == false) {
      baby.route.push_back(population[mother].route[j]);
      baby.cities_colected.push_back(population[mother].cities_colected[j]);
    }
  } // segunda metade da rota e coletas da mae
  return baby;
  // keep like that for now. test later heuristic for collecting bonus. search
  // prize collecting tsp for other possibilities
}// retorna um filho que terá metade inicial da rota do pai, segunda metade da
// rota da mãe (exceto cidades repetidas), coletando bônus caso eles coletem
// nessa cidade. NÃO CHAMA heuristica de carregamento de passageiros

Solution one_crossover(Instance instance, vector<Solution> &population) {
  int father = rand() % population.size(), mother = rand() % population.size();
  Solution baby = one_crossover(instance, population, father, mother);
  return baby;
}

vector<Solution> crossover(Instance instance, Population &parents) {
  vector<Solution> children;
  while (children.size() < parents.population.size()) {
    Solution baby = one_crossover(instance, parents.population);
    update_objectives(instance, baby);
    children.push_back(baby);
  }
  return children;
}
// retorna uma população de filhos igual ao tamanho da população de parentes

void save_data_routine(Generations& generations, BoundedParetoSet *EP, std::chrono::high_resolution_clock::time_point &inicio, int &biggest_multiple, int valuations, duration<double>& tempoTotal){
    auto fim = high_resolution_clock::now();
    duration<double> tempotrecho = duration_cast<duration<double>>(fim - inicio);//update tempo do trecho
    generations.durations.push_back(tempotrecho);
    tempoTotal += duration_cast<duration<double>>(fim - inicio); //somando tempo total
    Pareto_objectives sets_of_objectives = get_pareto_objectives(EP);
    generations.generations.push_back(sets_of_objectives);
    inicio = high_resolution_clock::now();
}

vector<Solution> crossover_and_mutate(Instance instance, Population &parents,int crossover_chance, int mutation_chance) {
  vector<Solution> children;
  while (children.size() < parents.population.size()) {
    Solution baby;
    if (crossover_chance < random_chance()) {
      baby = one_crossover(instance, parents.population);
      if (mutation_chance < random_chance()) {
        mutate_routes(instance, baby, 3);
        mutate_bonuses(instance, baby);
      }
      able_passengers(instance, baby);
      update_objectives(instance, baby);
    } else {
      baby = get_random_solution(instance); //esse aqui automaticamente dá update_objectives
    }
    if(solution_validity(instance, baby)){
      children.push_back(baby);
    }
  }
  return children;
}

vector<Solution> crossover_and_mutate_val1(Generations& generations,Instance instance, BoundedParetoSet *EP, Population &parents, int max_valuations, int crossover_chance, int mutation_chance, int& valuations,std::chrono::high_resolution_clock::time_point &inicio, int &biggest_multiple, duration<double>& tempoTotal) {
  vector<Solution> children;
  while (children.size() < parents.population.size()) {
    Solution baby;
    if (crossover_chance < random_chance()) {
      baby = one_crossover(instance, parents.population);
      able_passengers(instance, baby);
      update_objectives(instance, baby);
      if(solution_validity(instance, baby)){
        if(EP->adicionarSol(&baby)){
          continue;
        }
      }
      valuations++;
      int multiple = biggest_multiple+1;
      if(multiple*10000<=valuations){
        save_data_routine(generations,EP, inicio, biggest_multiple, valuations, tempoTotal);
        if(valuations>=max_valuations){
          return children;
        }
        biggest_multiple++;
      }
      if (mutation_chance < random_chance()) {
        mutate_routes(instance, baby, 3);
        mutate_bonuses(instance, baby);
        able_passengers(instance, baby);
        update_objectives(instance, baby);
        if(solution_validity(instance, baby)){
          if(EP->adicionarSol(&baby)){
            continue;
          }
        }
        valuations++;
        int multiple = biggest_multiple+1;
        if(multiple*10000<=valuations){
          save_data_routine(generations,EP, inicio, biggest_multiple, valuations, tempoTotal);
          if(valuations>=max_valuations){
            return children;
          }
          biggest_multiple++;
        }
      }
    } else {
      baby = get_random_solution(instance); //aqui dentro tem update objectives
      if(solution_validity(instance, baby)){
        if(EP->adicionarSol(&baby)){
          continue;
        }
      }
      valuations++;
      int multiple = biggest_multiple+1;
      if(multiple*10000<=valuations){
        save_data_routine(generations,EP, inicio, biggest_multiple, valuations, tempoTotal);
        if(valuations>=max_valuations){
          return children;
        }
        biggest_multiple++;
      }
    }
    if(solution_validity(instance, baby)){
      children.push_back(baby);
    }

  }
  return children;
}

void update_fronts(Population &population) {
  vector<vector<int>> solutions_dominated_by(population.population.size());
  vector<int> dominates(population.population.size(),0);
  vector<vector<int>> temp_fronts(population.population.size());
  int current_front = 0;

  for (int i = 0; i < population.population.size(); i++) {
    for (int j = 0; j < population.population.size(); j++) {
      if (i != j) {
        if (x_dominates_y(population.population[j], population.population[i])) {
          dominates[i]++;
          solutions_dominated_by[j].push_back(i);
        }
      }
    }
  }
   // Acima define quem domina quem e quantos dominam dada solução

  for (int i = 0; i < population.population.size(); i++) {
    if (dominates[i] == 0) {
      temp_fronts[0].push_back(i);
      dominates[i]=-1;
    }
  } //Define primeira fronte só com não-dominados

  while (!temp_fronts[current_front].empty()) {
    vector<int> next_front;
    for (int current_solution = 0; current_solution < temp_fronts[current_front].size();current_solution++) {//para cada solução x da fronte
      for (int j=0; j< solutions_dominated_by[temp_fronts[current_front][current_solution]].size();j++){ //para cada um dominado pela solução x
        //o numero da solução atual da current fronte: temp_fronts[current_front][current_solution]
        //soluções dominadas pela solução atual: solutions_dominated_by[temp_fronts[current_front][current_solution]] 
        dominates[solutions_dominated_by[temp_fronts[current_front][current_solution]][j]]--;
        if (dominates[solutions_dominated_by[temp_fronts[current_front][current_solution]][j]] == 0) {
          next_front.push_back(solutions_dominated_by[temp_fronts[current_front][current_solution]][j]);
        }
      }
    }
    current_front++;
    temp_fronts[current_front] = next_front;
  }
  population.fronts = temp_fronts;
  int i = 0;
  while (i < population.fronts.size()) {
    if (population.fronts[i].empty()) {
      population.fronts.erase(population.fronts.begin() + i);
    } 
    else{
    i++;
    }
  }//Limpando as frontes que ficaram vazias
}
// limpa as frontes atuais e as atualiza

void crowding_distance_assignment(Population &population) {

  int last = population.population.size();
  vector<pair<double, int>> objective_and_index(population.population.size());
  vector<double> crowding_distance(population.population.size());

  for (int objetivo = 1; objetivo <= 3; objetivo++) {
      if (objetivo == 1) { // cost
      for (int i = 0; i < last; i++) {
        objective_and_index[i].first = population.population[i].cost;
        objective_and_index[i].second = i;
      }//coloca os objetivos e os indices no vector de objetivos e indices
    } 
    else if (objetivo == 2) { // time
      for (int i = 0; i < last; i++) {
        objective_and_index[i].first = population.population[objective_and_index[i].second].time;
        objective_and_index[i].second = i;//isso ta certo?
      }
    }        
    else if (objetivo == 3) { // total bonus
      for (int i = 0; i < last; i++) {
        objective_and_index[i].first = population.population[objective_and_index[i].second].total_bonus;
        objective_and_index[i].second = i; // ta certo?
      }
    }
    heapSort(objective_and_index); //ordena crescente
    crowding_distance[objective_and_index[0].second] = 2147483647; // extremes receive biggest int possible
    crowding_distance[objective_and_index[last - 1].second] = 2147483647;        
    for (int i = 1; i < last - 1; i++) { // assigning crowding distance
      if(crowding_distance[objective_and_index[i].second]!=2147483647){
        crowding_distance[objective_and_index[i].second] += (objective_and_index[i + 1].first - objective_and_index[i - 1].first) /  (objective_and_index[last - 1].first - objective_and_index[0].first);
      }
    }
  }
  population.crowding_distance = crowding_distance;

}

int crowded_comparison_operator(Population &population, int solution1, int solution2) {
  int front_solution1, front_solution2;
  for (int front = 0; front < population.fronts.size(); front++) {
    for (int j = 0; j < population.fronts[front].size(); j++) {
      if (population.fronts[front][j] == solution1) {
        front_solution1 = front;
      }
      if (population.fronts[front][j] == solution2) {
        front_solution2 = front;
      }
    }
  }
  if (front_solution1 < front_solution2) {
    return solution1;
  } else if (front_solution1 > front_solution2) {
    return solution2;
  } else {
    if (population.crowding_distance[solution1] > population.crowding_distance[solution2]) {
      return solution1;
    } else {
      return solution2;
    }
  }
}

void kill_overpopulation_NSGA2(Population &population, int pop_max) { 
  int current_front = 0;
  vector<int> kill_list;
  vector<pair<double, int>> sorted_crowding_distance;
  while (current_front < population.fronts.size() and pop_max >= 0) {
    pop_max -= population.fronts[current_front].size();
    current_front++;
  } // após essa linha o current_front é o front em que começa a matança
  current_front--; // só pra voltar pra fronte correta
  for (int index : population.fronts[current_front]) {
    sorted_crowding_distance.push_back( make_pair(population.crowding_distance[index], index));
  } // após essa linha sorted_crowding_distance é uma lista com os crowding distances de cada um dos individuos do front atual e seus respectivos indices na fronte
  heapSort(sorted_crowding_distance);//heapsort em ordem crescente!!!

  for (int i = 0; i < (pop_max * -1) ; i++) { // erro aqui?
    kill_list.push_back(sorted_crowding_distance[i].second);
  } // após essa linha a kill_list tem as soluções que devem morrer da primeira fronte  
  current_front++;
  while (current_front < population.fronts.size()) {
    for (int i = 0; i < (population.fronts[current_front].size()); i++) {
      kill_list.push_back(population.fronts[current_front][i]);
    }
    current_front++;
  } // após essa linha já adicionou os individuos das frontes restantes na kill_list
  heapSort(kill_list); // heapsort ordena crescente
  for (int i = 0; i < kill_list.size(); i++){
    population.population.erase(population.population.begin() + kill_list[i]);
  }
}


Population NSGA2(Generations& NSGA2_generations, BoundedParetoSet *EP, duration<double>& tempoTotal, Instance instance, Population &population, int max_population, int max_valuations, int crossover_chance, int mutation_chance, int valuations_before_alg ) { // search tsp with profits
  // search heuristica for passageiros, bonus
  int valuations = valuations_before_alg ;
  int biggest_multiple = 0;
  Pareto_objectives sets_of_objectives = get_pareto_objectives(EP);
  NSGA2_generations.generations.push_back(sets_of_objectives);
  auto inicio = high_resolution_clock::now();
  while (valuations < max_valuations) {

    vector<Solution> children = crossover_and_mutate_val1(NSGA2_generations, instance, EP, population, max_valuations, crossover_chance, mutation_chance, valuations, inicio, biggest_multiple, tempoTotal);
    if(valuations>=max_valuations){
      break;
    }
    population.population.insert(population.population.end(), children.begin(), children.end());
    update_fronts(population);
    crowding_distance_assignment(population);
    kill_overpopulation_NSGA2(population, max_population);
  }
  population = get_non_dominated_population(population);
  return population;
}


Solution local_search_passengers(Instance instance, Solution &solution) {
  double initial_cost = solution.cost;

  for (int i = 0; i < instance.number_of_passengers - 1; i++) {
    for (int j = i + 1; j < instance.number_of_passengers; j++) {
      swap(solution.passengers_riding[i], solution.passengers_riding[j]);
      solution.cost = get_cost(instance, solution);

      if (check_passengers_riding(instance, solution) &&
          solution.cost < initial_cost) {
        return solution;
      } else {
        swap(solution.passengers_riding[i], solution.passengers_riding[j]);
        solution.cost = initial_cost;
        // if a better configuration is not found, swap back to initial
        // configuration
      }
    }
  }
  return solution;
}
// swaps passengers status. if a valid improvement is found or after testing all
// swaps, function stops

float calculate_euclidian_distance(vector<double> &weight_vector1, vector<double> &weight_vector2) {
  float distance = 0.0;
  for (int i = 0; i < weight_vector1.size(); i++) {
    distance += pow(weight_vector1[i] - weight_vector2[i], 2);
  }
  return sqrt(distance);
} // calcula distância entre dois vetores de pesos

vector<vector<float>> get_T_neighbors(vector<vector<double> > weight_vectors, int T_neighborhood_size){
    vector<vector<float>> T_neighbors(weight_vectors.size());
    for(int i = 0; i<weight_vectors.size();i++){
      vector<pair<double,int>> distances;
      for(int j = 0; j< weight_vectors.size();j++){
        if(i!=j){
          distances.push_back(make_pair(calculate_euclidian_distance(weight_vectors[i], weight_vectors[j]),j));
        }
      }
      heapSort(distances);
      for(int t = 0; t<T_neighborhood_size;t++){
        T_neighbors[i].push_back(distances[t].second);
      }
    }
    return T_neighbors;
}

vector<float> get_initial_z(Instance instance, Population population) {
  vector<float> z(3, 0);

  z[0] = population.population[0].cost;
  for (int i = 1; i < population.population.size(); i++) {
    double cost = get_cost(instance, population.population[i]);
    if (cost < z[0]) {
      z[0] = cost;
    }
  }

  z[1] = population.population[0].time;
  for (int i = 1; i < population.population.size(); i++) {
    double time = get_time(instance, population.population[i]);
    if (time < z[1]) {
      z[1] = time;
    }
  }

  z[2] = population.population[0].total_bonus;
  for (int i = 1; i < population.population.size(); i++) {
    double bonus = get_bonus(instance, population.population[i]);
    if (bonus > z[2]) {
      z[2] = bonus;
    }
  }

  return z;
} // retorna z = [cost, time, bonus]

void update_z(vector<float> z, Solution solution){
  if(solution.cost<z[0]){//step 2.3
    z[0] = solution.cost;
  }
  if(solution.time<z[1]){
    z[1] = solution.time;
  }
  if(solution.total_bonus> z[2]){
    z[2] = solution.total_bonus;
  }
}//z[cost,time,bonus]

double get_tchebycheff(Solution solution, vector<float> z, vector<double> weight_vector){
  double result = fabs(weight_vector[0]*(solution.cost - z[0]));
  if(result<fabs(weight_vector[1]*(solution.time - z[1]))){
    result = fabs(weight_vector[1]*(solution.time - z[1]));
  }
  if(result<fabs(weight_vector[2]*(solution.total_bonus - z[2]))){
    result = fabs(weight_vector[2]*(solution.total_bonus - z[2]));
  }
  return result;
}

pair<vector<Solution>, vector<int>> initial_population_moead (Population Population, vector<float> z, vector<vector<double> > weight_vectors){
  pair<vector<Solution>, vector<int>> initial_population;
  for(int i = 0; i<weight_vectors.size();i++){
    pair<double,int> best_for_this_vector;
    best_for_this_vector.second = 0;
    best_for_this_vector.first = get_tchebycheff(Population.population[0], z, weight_vectors[i]);

    for(int j=1;j<Population.population.size();j++){
      double j_tchebycheff = get_tchebycheff(Population.population[j], z, weight_vectors[i]);
      if(best_for_this_vector.first> j_tchebycheff){
        best_for_this_vector.second = j;
        best_for_this_vector.first = j_tchebycheff;
      }
      //definir qual da população é melhor para o vetor i
    }
    initial_population.first.push_back(Population.population[best_for_this_vector.second]);
    initial_population.second.push_back(best_for_this_vector.first);
  }
  return initial_population;
}

void update_neighboring_solutions(int n, pair<vector<Solution>, vector<int> > &internal_population, Solution baby, vector<vector<double> > weight_vectors, vector<vector<float>> T_neighbors,vector<float> z ){
  for(int t=0;t<T_neighbors[n].size(); t++){
    double baby_tchebycheff = get_tchebycheff(baby, z, weight_vectors[T_neighbors[n][t]]);
    double neighbor_tchebycheff = get_tchebycheff(internal_population.first[T_neighbors[n][t]], z, weight_vectors[T_neighbors[n][t]]);
    if(baby_tchebycheff < neighbor_tchebycheff){
        internal_population.first[T_neighbors[n][t]] = baby;
        internal_population.second[T_neighbors[n][t]] = baby_tchebycheff;
    }
  }
}



vector<vector<double>> generateWeightVectors(int NUMSUMPROBLEMAS) {
  vector<vector<double>> lambda(NUMSUMPROBLEMAS, vector<double>(3, 0.0f));
  int  s= 23; // NAO MUDAR
	int cont =0;
	double l1, l2, l3;
    for (int i=0; i<=s && cont<NUMSUMPROBLEMAS; i++){
        l1 = i;
        for (int j=0; j<=s-i && cont<NUMSUMPROBLEMAS; j++){
            l2 = j;
            l3 = s-i - l2;
            lambda[cont][0] = static_cast<double>(l1/s); // primeira coordenada
            lambda[cont][1] = static_cast<double>(l2/s); // segunda coordenada
            lambda[cont][2] = static_cast<double>(l3/s); // terceira coordenada
           
            // Esse é o seu vetor com três coordenadas normalizadas (a soma delas dá 1);
            //std::cout<<lambda[cont][0]<<" "<<lambda[cont][1]<<" " <<lambda[cont][2]<<std::endl;

            cont++; 
        }
    }
    return lambda;
}


Population MOEAD(Generations& MOEAD_generations, BoundedParetoSet *EP ,duration<double>& tempoTotal, Instance instance,  Population &population, int max_valuations, int N_subproblems, int T_neighborhood_size, vector<vector<double> > weight_vectors, int mutation_chance, int crossover_chance, int valuations_before_alg ) {
  //Population EP; //1.1
  vector<vector<float>> T_neighbors(weight_vectors.size());
  T_neighbors = get_T_neighbors(weight_vectors, T_neighborhood_size);//1.2
  vector<float> z = get_initial_z(instance, population);//1.4
  pair<vector<Solution>, vector<int>> internal_population = initial_population_moead(population, z, weight_vectors); //1.3
  //fim da inicialização
  int valuations = valuations_before_alg ;
  int biggest_multiple = 0;
  int multiple = biggest_multiple+1;
  Pareto_objectives sets_of_objectives = get_pareto_objectives(EP);
  MOEAD_generations.generations.push_back(sets_of_objectives);
  auto inicio = high_resolution_clock::now();
  while (valuations < max_valuations) {
    int n=0;
    while(n<N_subproblems){
        int father,mother;
        Solution baby;
        if(crossover_chance < random_chance()){
          father= T_neighbors[n][rand()%T_neighborhood_size];
          mother= T_neighbors[n][rand()%T_neighborhood_size];
          baby = one_crossover(instance, internal_population.first, father, mother);//2.1
          able_passengers(instance, baby);
          update_objectives(instance, baby);
          if(solution_validity(instance, baby)){
            if(EP->adicionarSol(&baby)){
              continue;
            }
          }
          valuations++;
          int multiple = biggest_multiple+1;
          if(multiple*10000<=valuations){
              save_data_routine(MOEAD_generations,EP, inicio, biggest_multiple, valuations, tempoTotal);
              if(valuations>=max_valuations){
                break;
              }
              biggest_multiple++;
          }
          
          if(mutation_chance < random_chance()){
            mutate_routes(instance,baby,3);
            mutate_bonuses(instance, baby);
            able_passengers(instance, baby);
            update_objectives(instance, baby);
            if(solution_validity(instance, baby)){
              if(EP->adicionarSol(&baby)){
                continue;
              }
            }
            valuations++;
            int multiple = biggest_multiple+1;
            if(multiple*10000<=valuations){
                save_data_routine(MOEAD_generations,EP, inicio, biggest_multiple, valuations, tempoTotal);
                if(valuations>=max_valuations){
                  break;
                }
                biggest_multiple++;
            }
          }//2.2
        }
        else{
          baby = get_random_solution(instance);
          if(solution_validity(instance, baby)){
            if(EP->adicionarSol(&baby)){
              continue;
            }
          }
          valuations++;
          int multiple = biggest_multiple+1;
          if(multiple*10000<=valuations){
              save_data_routine(MOEAD_generations,EP, inicio, biggest_multiple, valuations, tempoTotal);
              if(valuations>=max_valuations){
                break;
              }
              biggest_multiple++;
          }
        }
        if(solution_validity(instance,baby)==true){
          update_z(z,baby);//2.3
          update_neighboring_solutions(n,internal_population, baby, weight_vectors, T_neighbors, z); //2.4 
          n++;
        }
        //else{ cout<<"\n SOLUÇÂO ERRADA no MOEAD no fim do loop\n";}

    }    
  }
  return population;
}

double euclidian_distance_between_solutions(Solution solution1, Solution solution2) {
  double distance=0;
  distance += pow(solution1.cost - solution2.cost, 2);
  distance += pow(solution1.time - solution2.time, 2);
  distance += pow(solution1.total_bonus - solution2.total_bonus, 2);
  return sqrt(distance);
}

vector<int> solutions_strenght(vector<Solution> population){
  vector<int> strenght_vector(population.size());
  for(int i = 0; i<population.size();i++){
    int strenght = 0;
    for(int j = 0; j<population.size();j++){
      if(i!=j){
        if(x_dominates_y(population[i],population[j])){
          strenght++;
        }
      }
    }
    strenght_vector[i] = strenght;
  }
  return strenght_vector;
}

vector<int> raw_fitness_vector(vector<Solution> solutions, vector<int> strenght_vector){
  vector<int> raw_fitness_vector(solutions.size());
  strenght_vector = solutions_strenght(solutions);
  for(int i = 0; i<solutions.size();i++){
    for(int j = 0; j<solutions.size();j++){
      if(i!=j){
        if(x_dominates_y(solutions[j], solutions[i])){
          raw_fitness_vector[i] += strenght_vector[j];
        }
      }
    }
  }
  return raw_fitness_vector;
}

vector<vector<double>> distances_vector(vector<Solution> population){
  vector<vector<double> > distances_to_solution(population.size());
  for(int solution=0; solution<population.size(); solution++){
    vector<double> distances_to_neighbors(population.size());
    for(int neighbor = 0; neighbor<population.size(); neighbor++){
      if(solution!=neighbor){
        distances_to_neighbors[neighbor] = euclidian_distance_between_solutions(population[solution], population[neighbor]);
      }
    }
    heapSort(distances_to_neighbors);   
    distances_to_solution[solution] = distances_to_neighbors;
  }
  return distances_to_solution;
}

vector<double> density_vector(vector<Solution> population, vector<vector<double>> distances_vector){
  vector<double> density_vector(population.size());
  int k = sqrt(population.size());
  for(int solution = 0; solution<population.size(); solution++){
    density_vector[solution] = 1/(distances_vector[solution][k]+2);    
  }
  return density_vector;
}

vector<double> fitness_vector(vector<Solution> population, vector<int> raw_fitness, vector<vector<double>> distances_vector){
  vector<double> fitness_vector(population.size());
  vector<double> densities = density_vector(population, distances_vector);
  for(int i = 0; i<population.size();i++){
    fitness_vector[i] = raw_fitness[i]+ densities[i];
  }
  return fitness_vector;
}

Population SPEA2(Generations& SPEA2_generations, BoundedParetoSet *EP ,duration<double>& tempoTotal, Instance instance, Population population, int N_population_size, int N_archive_size, int max_valuations, int crossover_chance, int mutation_chance, int valuations_before_alg){
  Population A, archive, P = population, archive_and_P;
  Pareto_objectives sets_of_objectives = get_pareto_objectives(EP);
  SPEA2_generations.generations.push_back(sets_of_objectives);
  int valuations = valuations_before_alg ;
  int biggest_multiple = 0;
  auto inicio = high_resolution_clock::now();
  while(true){
    archive_and_P.population = P.population;
    archive_and_P.population.insert(archive_and_P.population.end(), archive.population.begin(), archive.population.end()); //adiciona a população do archive ao final da população
    archive_and_P.strenghts = solutions_strenght(archive_and_P.population);
    archive_and_P.raw_fitness = raw_fitness_vector(archive_and_P.population, archive_and_P.strenghts);
    archive_and_P.distances_to_other_solutions = distances_vector(archive_and_P.population);
    archive_and_P.fitness = fitness_vector(archive_and_P.population, archive_and_P.raw_fitness, archive_and_P.distances_to_other_solutions);

    pair<vector<Solution>,vector<int>> non_dominated_solutions;
    for(int solution=0; solution< archive_and_P.population.size(); solution++){
      if(archive_and_P.fitness[solution] <1){
        non_dominated_solutions.first.push_back(archive_and_P.population[solution]);
        non_dominated_solutions.second.push_back(solution);
      }
    }
    if(non_dominated_solutions.first.size() > N_archive_size){
      int excess = non_dominated_solutions.first.size() - N_archive_size;
      for(int solution=0; solution<excess; solution++){
        int solution_to_exclude=0;//0 faz diferença?
        for(int solution_to_compare=1; solution_to_compare<non_dominated_solutions.first.size(); solution_to_compare++){
          int k =0;
          while(archive_and_P.distances_to_other_solutions[non_dominated_solutions.second[solution_to_compare]][k] == archive_and_P.distances_to_other_solutions[non_dominated_solutions.second[solution_to_exclude]][k]){
            k++;
          }
          if(archive_and_P.distances_to_other_solutions[non_dominated_solutions.second[solution_to_compare]][k] < archive_and_P.distances_to_other_solutions[non_dominated_solutions.second[solution_to_exclude]][k]){
              solution_to_exclude = solution_to_compare;
          }
        }
        non_dominated_solutions.first.erase(non_dominated_solutions.first.begin()+solution_to_exclude);
        non_dominated_solutions.second.erase(non_dominated_solutions.second.begin()+solution_to_exclude);
        archive_and_P.distances_to_other_solutions = distances_vector(archive_and_P.population);
        //truncation technique tá feito, mas talvez com complexidade maior do que o necessário
      }
      archive.population = non_dominated_solutions.first;
    }
    else if(non_dominated_solutions.first.size() < N_archive_size){
      int missing = N_archive_size - non_dominated_solutions.first.size();
      archive.population = non_dominated_solutions.first;
      vector<pair<double, int>> dominated_solutions;
      for(int solution=0; solution< archive_and_P.population.size(); solution++){
        if(archive_and_P.fitness[solution] >=1){
          dominated_solutions.push_back(make_pair(archive_and_P.fitness[solution], solution));
        }
      }
      heapSort(dominated_solutions);
      for(int inclusion=0; inclusion<missing; inclusion++){        
        archive.population.push_back(archive_and_P.population[dominated_solutions[inclusion].second]);
      }
    }
    else{
      archive.population = non_dominated_solutions.first;
    }
    if(valuations>=max_valuations){
      A.population = non_dominated_solutions.first;
      break;
    }
    else{
      P.population = crossover_and_mutate_val1(SPEA2_generations, instance, EP, population, max_valuations, crossover_chance, mutation_chance, valuations, inicio, biggest_multiple, tempoTotal);
    }//checar aqui de novo
  }

  return A;
}

void MO_vns(Instance instance, Population &population) {
  int k_max = 4;
  for (int i = 0; i < population.population.size(); ++i) {
    Solution current_solution = population.population[i];
    Solution new_solution = current_solution;
    for (int k = 1; k <= k_max; ++k) {
      switch (k) {
      case 1:
        mutate_routes(instance, new_solution, 0);
        able_passengers(instance, new_solution);
        break;
      case 2:
        mutate_routes(instance, new_solution, 1);
        able_passengers(instance, new_solution);
        break;
      case 3:
        mutate_routes(instance, new_solution, 2);
        able_passengers(instance, new_solution);
        break;
      case 4:
        invert_bonuses(instance, new_solution);
        able_passengers(instance, new_solution);
        break;
      }
      if (x_dominates_y(new_solution,current_solution)) { // Verificar se a nova solução
                                             // domina a solução original
        population.population[i] = new_solution;
        break; // Sair do loop k
      }
    }
  }
}

ObjectiveBounds computeObjectiveBounds(const vector<Solution>& population) {
  ObjectiveBounds bounds;
  
  bounds.min_cost = bounds.min_time = numeric_limits<double>::max();
  bounds.max_cost = bounds.max_time = numeric_limits<double>::lowest();
  bounds.min_bonus = bounds.max_bonus = 0;  // O bônus nunca será negativo

  for (const auto& sol : population) {
      if (sol.cost < bounds.min_cost) bounds.min_cost = sol.cost;
      if (sol.cost > bounds.max_cost) bounds.max_cost = sol.cost;
      
      if (sol.time < bounds.min_time) bounds.min_time = sol.time;
      if (sol.time > bounds.max_time) bounds.max_time = sol.time;
      
      if (sol.total_bonus > bounds.max_bonus) bounds.max_bonus = sol.total_bonus;
  }

  return bounds;
}

NormalizedObjectives scaleObjectives(const Solution& sol, const ObjectiveBounds& bounds, double rho) {
  NormalizedObjectives norm;

  // Ajuste do limite superior
  double adjusted_max_cost = bounds.max_cost + rho * (bounds.max_cost - bounds.min_cost);
  double adjusted_max_time = bounds.max_time + rho * (bounds.max_time - bounds.min_time);

  // Normalização no intervalo [0,1]
  norm.cost = (sol.cost - bounds.min_cost) / (adjusted_max_cost - bounds.min_cost);
  norm.time = (sol.time - bounds.min_time) / (adjusted_max_time - bounds.min_time);
  norm.bonus = (sol.total_bonus - bounds.min_bonus) / (bounds.max_bonus - bounds.min_bonus + 1e-9); // Evita divisão por zero

  return norm;
}

vector<NormalizedObjectives> normalizePopulation(const vector<Solution>& population, double rho) {
  ObjectiveBounds bounds = computeObjectiveBounds(population);
  vector<NormalizedObjectives> normalizedPopulation;
  normalizedPopulation.reserve(population.size());

  for (const auto& sol : population) {
      normalizedPopulation.push_back(scaleObjectives(sol, bounds, rho));
  }

  return normalizedPopulation;
}


double calculateIndicator(const NormalizedObjectives& x1, const NormalizedObjectives& x2) {
  return (x2.cost - x1.cost) + (x2.time - x1.time) - (x2.bonus - x1.bonus);
}

vector<vector<double>> computeIndicatorMatrix(const vector<NormalizedObjectives>& normalized_population) {
  size_t pop_size = normalized_population.size();
  vector<vector<double>> indicator_matrix(pop_size, vector<double>(pop_size, 0.0));

  for (size_t i = 0; i < pop_size; ++i) {
      for (size_t j = 0; j < pop_size; ++j) {
          if (i != j) {
              indicator_matrix[i][j] = calculateIndicator(normalized_population[i], normalized_population[j]);
          }
      }
  }

  return indicator_matrix;
}

double computeMaxIndicatorValue(const vector<vector<double>>& indicator_matrix) {
  double maxAbsIndicator = 0.0;
  for (size_t i = 0; i < indicator_matrix.size(); ++i) {
      for (size_t j = 0; j < indicator_matrix[i].size(); ++j) {
          if (i != j) {
              maxAbsIndicator = std::max(maxAbsIndicator, std::abs(indicator_matrix[i][j]));
          }
      }
  }
  return maxAbsIndicator;
}

vector<double> computeFitness(const vector<NormalizedObjectives>& normalized_population, double kappa) {
  size_t pop_size = normalized_population.size();
  vector<double> fitness_values(pop_size, 0.0);
  
  // Calcula a matriz de indicadores
  vector<vector<double>> indicator_matrix = computeIndicatorMatrix(normalized_population);
  
  // Obtém o maior valor absoluto do indicador
  double max_indicator = computeMaxIndicatorValue(indicator_matrix);

  // Computa o fitness de cada indivíduo
  for (size_t i = 0; i < pop_size; ++i) {
      double fitness = 0.0;
      for (size_t j = 0; j < pop_size; ++j) {
          if (i != j) {
              fitness += -exp(-indicator_matrix[j][i] / (max_indicator * kappa));
          }
      }
      fitness_values[i] = fitness;
  }

  return fitness_values;
}
int findWorstIndividual(const vector<double>& fitness_values) {
  int worst_index = -1;
  double worst_fitness = numeric_limits<double>::max();

  for (size_t i = 0; i < fitness_values.size(); ++i) {
      if (fitness_values[i] < worst_fitness) {
          worst_fitness = fitness_values[i];
          worst_index = i;
      }
  }
  return worst_index;
}

void environmentalSelection(vector<Solution>& population, vector<double>& fitness_values, int MU) {
  while (population.size() > MU) {
      int worst_index = findWorstIndividual(fitness_values);
      if (worst_index != -1) {
          // Remove o pior indivíduo da população e do vetor de fitness
          population.erase(population.begin() + worst_index);
          fitness_values.erase(fitness_values.begin() + worst_index);
      }
  }
}


Population IBEA(Generations& IBEA_generations, BoundedParetoSet *EP,Instance instance, Population received_population, int N_population_size, int max_valuations,  int mutation_chance, int crossover_chance, int valuations_before_alg, int MU, double kappa, double scalingFactor, double rho){
    Pareto_objectives sets_of_objectives = get_pareto_objectives(EP);
    IBEA_generations.generations.push_back(sets_of_objectives);
    Population population = received_population;
    int valuations = valuations_before_alg;
    int biggest_multiple = 0;

    while (valuations < max_valuations) {
        vector<NormalizedObjectives> normalized_population = normalizePopulation(population.population, rho);
        vector<double> fitness_values = computeFitness(normalized_population,kappa);
        environmentalSelection(population.population, fitness_values, MU);
        
        vector<Solution> P;
        P = crossover_and_mutate(IBEA_generations, instance, EP, population, max_valuations, crossover_chance, mutation_chance, valuations, biggest_multiple);
        population.population.insert(population.population.end(), P.begin(), P.end()); 
        if(valuations > max_valuations){
          break; 
        }
        cout<<"IBEA terminou crossover and mutate com "<<P.size()<<" filhos, valuations:"<<valuations<<endl;
    }
    population = get_non_dominated_population(population);
    return population;
}

void print_solution_NSGA2(Instance instance, Solution &solution) {
  cout << "Rota: ";
  for(int i = 0; i < solution.route.size(); i++){
    cout << solution.route[i] << " ";
  }
  cout<<endl;

  cout<<"custos da rota:"<<endl;
  for(int city = 0; city < solution.route.size(); city++){
    int next_city;
    if(city == solution.route.size() - 1){
      next_city = 0;
    }
    else{
      next_city = city+1;
    }
    cout << instance.cost_matrix[solution.route[city]][solution.route[next_city]]<< " ";
  }
  cout<<endl;

  cout<<"Tempos da rota:"<<endl;
  for(int city = 0; city < solution.route.size(); city++){
    int next_city;
    if(city == solution.route.size() - 1){
      next_city = 0;
    }
    else{
      next_city = city+1;
    }
    cout << instance.time_matrix[solution.route[city]][solution.route[next_city]]<< " ";
  }
  cout<<endl;

  cout << "Cities_colected: "<<endl;
  bool any_city_colected = false;
  for(int i =0; i<solution.cities_colected.size();i++){
    if(solution.cities_colected[i]==true){
      cout<<"City: "<< solution.route[i] << " bonus " << instance.bonus_and_time[solution.route[i]].first << " tempo " << instance.bonus_and_time[solution.route[i]].second << endl;
      any_city_colected = true;
    }
  }
  if(!any_city_colected){
    cout<<"Nenhuma cidade coletada"<<endl;
  }
  cout<<endl;
}

void print_population(Instance instance, Population population) {
  for (int i = 0; i < population.population.size(); i++) {
    cout << "Solução " << i << ":" << endl;
    print_solution(instance, population.population[i]);
  }
  cout << "Quantidade de soluções: " << population.population.size() << endl;
}

void shorter_print_population(Instance instance, Population population){
  cout<<"RESUMIDO (solução: custo tempo bonus)\n";
  for(int i =0;i<population.population.size();i++){
    cout<<population.population[i].cost<<" "<<population.population[i].time<<" "<<population.population[i].total_bonus<<" "<<endl;
  }
  cout<<"tamanho da população: "<<population.population.size()<<endl;
}

void print_instance(Instance instance){
  cout << "Number of cities: " << instance.number_of_cities << endl;
  cout << "Number of passengers: " << instance.number_of_passengers << endl;
  cout << "Car capacity: " << instance.car_capacity <<"\n " <<endl;
}

void compare_pop(vector<Solution> solutions1, vector<Solution> solutions2, bool printtests){
  int contador_de_pop_1_so_domina=0;
  int contador_de_pop_1_so_dominado=0;
  int contador_de_indiferentes=0;
  int contador_de_pop_1_sanduiche=0;

  for(int i=0;i<solutions1.size();i++){
    bool bandeira_domina = false;
    bool bandeira_dominado = false;
    for(int j =0;j<solutions2.size();j++){
      if(x_dominates_y(solutions1[i], solutions2[j])){
          bandeira_domina =true;
          if(printtests == true){
            cout<<"\n valores do vetor que domina: "<<solutions1[i].cost<<" "<<solutions1[i].time<<" "<<solutions1[i].total_bonus<<endl;
            cout<<"valores do vetor dominado: "<<solutions2[j].cost<<" "<<solutions2[j].time<<" "<<solutions2[j].total_bonus<<"\n"<<endl;
          }
      }
      else if(x_dominates_y(solutions2[j], solutions1[i])){
          bandeira_dominado =true;

      }
    }

    if(bandeira_domina ==true and bandeira_dominado==false){
      contador_de_pop_1_so_domina++;
    }
    else if(bandeira_domina ==false and bandeira_dominado ==true ){
      contador_de_pop_1_so_dominado++;
    }
    else if(bandeira_dominado == false and bandeira_domina ==false){
      contador_de_indiferentes++;
    }
    else{
      contador_de_pop_1_sanduiche++; //quando domina e é dominado
    }
  }
  cout<<"População 1 com tamanho "<<solutions1.size()<<" e população 2 com tamanho "<<solutions2.size()<<endl;
  cout<<"A população 1 teve "<<contador_de_pop_1_so_domina<<" soluções não-dominadas dominando ao menos um da população 2" <<endl;
  cout<<"A população 1 teve "<<contador_de_pop_1_so_dominado<<" soluções dominadas pela população 2 sem dominar nenhum" <<endl;
  cout<<"A população 1 teve "<<contador_de_indiferentes<<" soluções não-dominadas que não dominam" <<endl;  
  cout<<"A população 1 teve "<<contador_de_pop_1_sanduiche<<" soluções dominadas que também dominam" <<endl;

}

vector<Solution> matrix_to_solution(vector<vector<double>> matriz){
    vector<Solution> solutions;

    for (int i = 0; i < matriz.size(); ++i) {
        Solution sol;
        sol.cost = matriz[i][0];
        sol.time = static_cast<int>(matriz[i][1]);
        sol.total_bonus = static_cast<int>(matriz[i][2]);
        solutions.push_back(sol);
    }
    return solutions;
}

void print_weight_vector(vector<vector<float>> weight_vectors){
  for(int i =0; i<weight_vectors.size();i++){
    cout<<"vetor "<<i<<": ";
    for(int j =0; j<weight_vectors[i].size();j++){
      cout<<weight_vectors[i][j]<<"   ";
    }
    cout<< "\n";
  }  
}

void create_all_directories(int num_tests, int num_instances, std::string folder_name, vector<string> instances_addresses) {
  for (const std::string& algo : {"NSGA2", "MOEAD"}) { 
    std::string algo_dir = folder_name + "/" + algo;
    for (int instance = 0; instance < num_instances; ++instance) {
        std::string instance_dir = algo_dir + "/"+instances_addresses[instance];  
          for (int test = 0; test < num_tests; ++test) {
            std::string teste_dir = instance_dir + "/teste_" + std::to_string(test);
            filesystem::create_directories(teste_dir);
           }
        }
    }
}

void fill_complete_archive(Instance instance, ofstream &arquivo_completo, Generations saved_generations, int generation, int pareto_set){
    arquivo_completo<<"Rota: ";
    for(int i = 0; i < saved_generations.generations[generation].solutions[pareto_set].route.size(); i++){
      arquivo_completo<< saved_generations.generations[generation].solutions[pareto_set].route[i] << " ";
    }
    arquivo_completo<<endl;
    arquivo_completo<<"custos da rota: ";
    for(int city = 0; city < saved_generations.generations[generation].solutions[pareto_set].route.size(); city++){
    int next_city;
    if(city == saved_generations.generations[generation].solutions[pareto_set].route.size() - 1){
      next_city = 0;
    }
    else{
      next_city = city+1;
    }
    arquivo_completo << instance.cost_matrix[saved_generations.generations[generation].solutions[pareto_set].route[city]][saved_generations.generations[generation].solutions[pareto_set].route[next_city]]<< " ";
  }
  arquivo_completo<<endl;

  arquivo_completo<<"Tempos da rota: ";
  for(int city = 0; city < saved_generations.generations[generation].solutions[pareto_set].route.size(); city++){
    int next_city;
    if(city == saved_generations.generations[generation].solutions[pareto_set].route.size() - 1){
      next_city = 0;
    }
    else{
      next_city = city+1;
    }
    arquivo_completo<< instance.time_matrix[saved_generations.generations[generation].solutions[pareto_set].route[city]][saved_generations.generations[generation].solutions[pareto_set].route[next_city]]<< " ";
  }
  arquivo_completo<<endl;

  arquivo_completo<< "Cities_colected: "<<endl;
  bool any_city_colected = false;
  for(int i =0; i<saved_generations.generations[generation].solutions[pareto_set].cities_colected.size();i++){
    if(saved_generations.generations[generation].solutions[pareto_set].cities_colected[i]==true){
      arquivo_completo<<"City: "<< saved_generations.generations[generation].solutions[pareto_set].route[i] << " bonus " << instance.bonus_and_time[saved_generations.generations[generation].solutions[pareto_set].route[i]].first << " tempo " << instance.bonus_and_time[saved_generations.generations[generation].solutions[pareto_set].route[i]].second << endl;
      any_city_colected = true;
    }
  }
  if(!any_city_colected){
    arquivo_completo<<"Nenhuma cidade coletada"<<endl;
  }
  int passengers_on = 0;
  for (int i = 0; i < saved_generations.generations[generation].solutions[pareto_set].passengers_riding.size(); i++) {
    if (saved_generations.generations[generation].solutions[pareto_set].passengers_riding[i]) {
      passengers_on++;
      arquivo_completo<< "passageiro "<<i<<" com origem "<<instance.passengers[i].origin<<", destino: "<<instance.passengers[i].destiny<< " custo max: "<<instance.passengers[i].max_cost<<" tempo max: "<<instance.passengers[i].max_time<<endl;
    }
  }
  arquivo_completo << "Passageiros embarcados: " << passengers_on << endl;
  arquivo_completo<< "Custo: " << saved_generations.generations[generation].solutions[pareto_set].cost << endl;
  arquivo_completo<< "Bonus: " << saved_generations.generations[generation].solutions[pareto_set].total_bonus << endl;
  arquivo_completo<< "Tempo: " << saved_generations.generations[generation].solutions[pareto_set].time << "\n \n";

}

void create_test_archives(Instance instance, string algorithm_name, string folder_name, vector<string> instances_addresses, int current_instance, int current_test, Generations saved_generations, unsigned int seed, duration<double> tempoTotal){
  for(int generation = 0; generation < saved_generations.generations.size(); generation++){
    string endereço_do_arquivo = folder_name+ "/"+algorithm_name+"/" + instances_addresses[current_instance] +"/"+ "teste_" + to_string(current_test) +"/paretogeracao_" + to_string(generation) + ".txt";
    string endereço_do_arquivo_completo = folder_name+ "/"+algorithm_name+"/" + instances_addresses[current_instance] +"/"+ "teste_" + to_string(current_test) +"/paretogeracao_" + to_string(generation) + "_completo.txt";
    string endereço_dados = folder_name+ "/"+algorithm_name+"/" + instances_addresses[current_instance] +"/"+ "teste_" + to_string(current_test) + "/dados";
    ofstream arquivo(endereço_do_arquivo);
    ofstream arquivo_completo(endereço_do_arquivo_completo);
    ofstream arquivo_dados(endereço_dados);
    for(int pareto_set = 0; pareto_set< saved_generations.generations[generation].pareto_set.size() ; pareto_set++){
      arquivo<<saved_generations.generations[generation].pareto_set[pareto_set].cost<<" "<< saved_generations.generations[generation].pareto_set[pareto_set].time<<" " << saved_generations.generations[generation].pareto_set[pareto_set].total_bonus;
      fill_complete_archive(instance, arquivo_completo, saved_generations, generation, pareto_set);
      //adicionar rota, cities, passengers, custos, tempos, etc
      if(pareto_set<saved_generations.generations[generation].pareto_set.size()-1){
        arquivo<<endl;
        arquivo_completo<<endl;
      }
    }
    arquivo_dados<< "seed: "<<seed<<endl;
    arquivo_dados<< "tempo total em segundos: "<<tempoTotal.count()<<endl;
    for(int duration =0;duration < saved_generations.durations.size(); duration++){
      arquivo_dados<<"tempo da geração "<<duration<< ":"<<endl<< saved_generations.durations[duration].count()<<endl;
    }
    arquivo.close();
    arquivo_completo.close();
    arquivo_dados.close();
  }
}

void create_archive_IRACE(Pareto_objectives sets_of_objectives){
  string endereço_do_arquivo =  "testeirace.txt";
  ofstream arquivo(endereço_do_arquivo);
  for(int pareto_set = 0; pareto_set< sets_of_objectives.pareto_set.size() ; pareto_set++){
    arquivo<<sets_of_objectives.pareto_set[pareto_set].cost<<" "<< sets_of_objectives.pareto_set[pareto_set].time<<" " << sets_of_objectives.pareto_set[pareto_set].total_bonus;
    if(pareto_set<sets_of_objectives.pareto_set.size()-1){
      arquivo<<endl;
    }
  }
  arquivo.close();
}

int main(int argc, char* argv[]) {
  //
  //SETUP 
  int max_population = 300;
  int N_subproblems = 300;
  int max_valuations = 80000;
  Instance instance;
  int crossover_chance;
  int mutation_chance;

  //MOEAD
  int T_neighborhood_size ;
  vector<vector<double>> weight_vectors = generateWeightVectors(N_subproblems);

  //IBEA
  double rho = 0.1; // Parâmetro para ajuste do limite superior
  double scalingFactor = 0.1; // Fator de escala (ρ)
  double kappa = 0.05;       // Fator de ajuste para a aptidão
  int MU = 300;              // Tamanho final da população

  //int N_population_size = 100; 
  for(int i =0; i<argc; i++){
    string teste = argv[i];
    if(teste == "-i" ){
      instance = get_instance(argv[i+1]);
    }
    if(teste == "--crossover"){
      crossover_chance = stof(argv[i+1]);
    }
    if(teste == "--mutation"){
      mutation_chance = stof(argv[i+1]);
    }
    if(teste == "--scalingFactor"){
      scalingFactor = stof(argv[i+1]);
    }
    if(teste == "--kappa"){
      kappa = stof(argv[i+1]);
    }
    if(teste == "--rho"){
      rho = stof(argv[i+1]);
    }
  }
    //Testes IBEA
      BoundedParetoSet *EP = new BoundedParetoSet(); 
      Population population;
      Generations saved_generations;
      unsigned int seed = static_cast<unsigned int>(time(0)); // pra salvar a seed depois
      srand(seed);

      int valuations_before_alg =0;
      Population initial_population = get_random_population(instance, EP, max_population,valuations_before_alg); // inicializar população
      int max_valuations_after_alg = max_valuations- valuations_before_alg;
      duration<double> tempoTotal = duration<double>::zero();//Variável para armazenar o tempo total
      population = IBEA(saved_generations, EP, instance, initial_population, max_population, max_valuations_after_alg, mutation_chance, crossover_chance, valuations_before_alg, MU, kappa, scalingFactor, rho);
      
      Pareto_objectives sets_of_objectives = get_pareto_objectives(EP);
      delete EP;
      
      create_archive_IRACE(sets_of_objectives);
      std::ofstream arquivo("iracehyp.txt", std::ios::out);
      arquivo.close();
      string comando;
      comando = "hyp.exe hyp_ind_param_irace.txt testeirace.txt testeirace.txt iracehyp.txt";
      int result = system(comando.c_str());
      std::ifstream arquivoirace("iracehyp.txt");
      std::string linha;
      while (std::getline(arquivoirace, linha)) {
          //linha = "-"+ linha;
          double linha2 = stof(linha);
          cout <<linha2;
          //int numero = 24;
          //cout<<numero;
      }
      arquivoirace.close();
}

#endif