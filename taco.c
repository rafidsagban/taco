#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <sys/param.h>

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
#define p 0.1
#define N 100
#define EPS_REL 0.1
#define MAX_ITER 100
#define COUNT 101

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
int F(Ant* a, Ant* b, int f_idx);
int Compare(Ant* a, Ant* b, int* order, int len);
void simulate();
void get_best_ant();
void update_tau();

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
Ant* loc_best;

int main(int argc, char** argv){
    srand(time(NULL));
    read_input();
    init();
    create_feasible_solution();
    simulate();
    int i;
    for(i=0;i<loc_best->route_len;i++)
        printf("%d->", loc_best->route[i]);

    free_mem();
    return 0;
}

void read_input(){
    char buf[128];
    int i,j;
    scanf("%s\n", (char*)&name);
    
    for(i=0;i<3;i++)
        scanf("%s \n", buf);
    
    scanf("%d %lf", &count, &capacity);
    count = COUNT;
    
    for(i=0;i<12;i++)
        scanf("%s \n", buf);

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

void simulate(){
    int i=0;
    do{
        create_feasible_solution();
        get_best_ant();
        printf("%d routes, %6.3f delay, %6.3f overload, %6.3f length\n", 
                loc_best->route_count, loc_best->max_delay, 
                loc_best->max_overload, loc_best->total_time);
        update_tau();
    } while (++i<MAX_ITER);   
}

void get_best_ant(){
    int i;
    loc_best = &(ants[0]);
    int order[4] = {3,1,2,0};
    for(i=1;i<N;i++){
        if(Compare(&(ants[i]), loc_best, (int*)&order, 4) < 0)
            loc_best = &(ants[i]);
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

void update_tau(){
    int i,j, from, to;
    Ant* ant;
    for(i=0;i<count*count;i++)
        tau[i] = (1-p) * tau[i];
    
    ant = loc_best;
    for(j=0;j<ant->route_len-1;j++){
        from = ant->route[j];
        to = ant->route[j+1];
        tau[from*count+to] += 3000.0/ant->total_time;
    }
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

int Compare(Ant* a, Ant* b, int* order, int len){
    int res = 0;
    if(len>0){
        res = F(a,b,order[0]);
        if(res == 0 && len > 1)
            res = Compare(a,b,&(order[1]), len-1);
    }
    return res;
}

int F(Ant* a, Ant* b, int f_idx){
    double va = 0, vb = 0;
    switch(f_idx){
        case 0:
            va = a->route_count;
            vb = b->route_count;
            break;
        case 1:
            va = a->max_delay;
            vb = b->max_delay;
            break;
        case 2:
            va = a->max_overload;
            vb = b->max_overload;
            break;
        case 3:
            va = a->total_time;
            vb = b->total_time;
            break;
    }
    if(va-vb < -EPS_REL*MIN(va,vb))
        return -1;
    else if(va-vb > EPS_REL*MIN(va,vb))
        return 1;
    return 0;
}