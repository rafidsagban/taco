#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>

typedef struct _twnd {
    double start;
    double stop;
} TimeWindow;
typedef struct _coord {
    double x;
    double y;
} Coord;
typedef struct _ant {
	int route_count;
	double route_len;
	double total_time;
	double max_delay;
	double max_overload;
	int* route;
} Ant;

#define alpha 1.0
#define beta 1.0
#define N 10
#define EPS_REL 0.05

inline double distance(Coord* p1, Coord* p2);
inline double max_overload(Ant* ant);
void read_input();
void free_mem();
double my_rand();
int roulette_next(double* ar, int len);
void init();
void create_ants();
void delete_ants();
void ant_action(int idx);
void create_feasible_solution();

#if !DEBUG

Coord* coords;
TimeWindow* tw;
double* w;
double* svc;
double* d;
double* tau;
Ant* ants;
int* IN;

#else

#define COUNT 26
Coord coords[COUNT];
TimeWindow tw[COUNT];
double w[COUNT];
double svc[COUNT];
double d[COUNT*COUNT];
double tau[COUNT*COUNT];
Ant ants[N];
int IN[COUNT];

#endif

char name[10];
int count;
double capacity;

int main(int argc, char** argv){
    srand(time(NULL));
    read_input();
    init();
    create_feasible_solution();

    free_mem();
    return 0;
}

void read_input(){
    scanf("%s\n", (char*)&name);
    scanf("%d %lf", &count, &capacity);
    count++;

#if !DEBUG
    coords = (Coord*)malloc(sizeof(Coord)*count);
    w = (double*)malloc(sizeof(double)*count);
    svc = (double*)malloc(sizeof(double)*count);
    d = (double*)malloc(sizeof(double)*count*count);
    tau = (double*)malloc(sizeof(double)*count*count);
    tw = (TimeWindow*)malloc(sizeof(TimeWindow)*count);
    ants = (Ant*)malloc(sizeof(Ant)*N);
    IN = (int*)malloc(sizeof(int)*count);
#endif

    int i,j;
    for(i=0;i<count;i++){
        scanf("%d", &j);
        scanf("%lf %lf %lf", &(coords[j].x), &(coords[j].y), &(w[j]));
        scanf("%lf %lf %lf", &(tw[j].start), &(tw[j].stop), &(svc[j]));
    }

    for(i=0;i<count;i++)
        for(j=i+1;j<count;j++)
            d[i*count+j]=d[j*count+i]=distance(&(coords[i]),&(coords[j]));
}

void free_mem(){
	delete_ants();

#if !DEBUG
    free(coords);
    free(w);
    free(svc);
    free(d);
    free(tau);
    free(tw);
    free(ants);
    free(IN);
#endif
}

inline double distance(Coord* p1, Coord* p2){
    return sqrt((p2->x-p1->x)*(p2->x-p1->x)+(p2->y-p1->y)*(p2->y-p1->y));
}

inline double max_overload(Ant* ant){
	int* route = (ant->route);
	int i;
	double weight = 0, max = 0, ovl = 0;
	for(i=0;i<ant->route_len;i++){
		if(i>0 && route[i] == 0){
			ovl = weight - capacity;
			weight = 0;
			if(ovl > max)
				max = ovl;
		} else {
			weight += w[route[i]];
		}
	}
	return max;
}

int roulette_next(double* ar, int len){
    double s = my_rand();
    int i;
    for(i=0;i<len;s-=ar[i],i++)
        if(s<ar[i])
            return i;
    return -1;
}

double my_rand(){
    return ((double)rand())/RAND_MAX;
}

void init(){
	double tau0 = 1.0/count;
	int i;
	for(i=0;i<count*count;i++)
		tau[i] = tau0;

    create_ants();
}

void create_ants(){
	int i;
	for(i=0;i<N;i++){
		ants[i].route_count = 0;
		ants[i].route_len = 0;
		ants[i].total_time = 0;
		ants[i].route = (int*)malloc(sizeof(int)*(count*2 - 1));
	}
}

void delete_ants(){
	int i;
	for(i=0;i<N;i++){
		free(ants[i].route);
	}
}

void create_feasible_solution(){
	int i;
	for(i=0;i<N;i++)
		ant_action(i);
}

void ant_action(int idx){
	int i, cur = 0, route_len=0, flag;
	double time = 0, dist, a, b, c, sum, n,
			total_time = 0, route_count = 0,
			max_delay = 0;

	Ant* ant = &(ants[idx]);
	int* visited = (int*)malloc(sizeof(int)*count);
	int* route = (int*)malloc(sizeof(int)*count*2);
	double* P = (double*)malloc(sizeof(double)*count);
	memset(visited, 0, sizeof(int)*count);
	memset(P, 0, sizeof(double)*count);
	do{

		do{
			route[route_len++] = cur;
			visited[cur] = 1;
			sum = 0;
			for(i=0;i<count;i++){
				if(visited[i] || tw[i].stop < time)
					P[i] = 0;
				else {
					a = time + d[cur*count+i];
					b = a > tw[i].start ? a : tw[i].start;
					c = a > b ? a : b;
					dist = (tw[i].stop - time)*(c-time);
					n = 1 > dist - IN[i] ? 1 : dist - IN[i];
					n = 1.0/n;
					P[i] = pow(tau[cur*count+i], alpha)*pow(n, beta);
					sum += P[i];

					if(tw[i].start - a > max_delay)
						max_delay = tw[i].start - a;
				}
			}
			if(sum == 0) break;
			for(i=0;i<count;i++)
				P[i] = P[i] / sum;

			int sel = roulette_next(P, count);
			if(sel == -1) break;

			time += d[cur*count+sel] + svc[sel];
			cur = sel;

		} while(1);
		total_time += time;
		route_count++;

		cur = 0;
		time = 0;

		flag = 0;
		for(i=0;i<count;i++)
			if(!visited[i])
				flag = 1;

	} while(flag);

	route[route_len++] = 0;

	memcpy(ant->route, route, sizeof(int)*(2*count-1));
	ant->route_len = route_len;
	ant->route_count = route_count;
	ant->total_time = total_time;
	ant->max_delay = max_delay;
	ant->max_overload = max_overload(ant);

	free(P);
	free(visited);
	free(route);
}
