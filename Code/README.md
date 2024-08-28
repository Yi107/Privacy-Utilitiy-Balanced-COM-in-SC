# Code Description

Data_structure.h initializes the strucure of Worker and Request used in this code.

Data_structure.cpp defines some basic functions about these structures.

geo.cpp and geo.h define the basic functions to calculate the distance between two points based on their longitude and latitude. It also defines the functions which transfers the longitude and latitude to x-y coordinate axis. 

privacy.cpp defines all perturbation functions and random number generation functions.

Match.cpp defines all matching algorithms and pricing algorithms. The main matching algorithm is basedMatch(worker, tasks,matching_results,matching_number, revenue, pricing_epsilon,geo_epsilon).

Main.cpp defines the final running process. read_file(worker_filename, task_filename, workers_set, tasks_set,worker_radius).

Readers could perform int main() algorithm with different worker_file, task_file, worker_radius, pricing_epsilon,geo_epsilon as input to get different results.

