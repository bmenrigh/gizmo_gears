allocatemem();
allocatemem();
allocatemem();
allocatemem();
allocatemem();

setrand(extern("date +%s"));

/* param_n = 12; */
/* param_r = sqrt(2); */
/* param_r = 1.377; */

/* param_n = 7; */
/* param_r = 1.623579; */
/* param_r = sqrt(2.63600876924); */
/* param_r = sqrt(3); */
/* param_r = 1 + (2.0 / 3.0); */
/* param_r = 1.6275; */

/* param_n = 8; */
/* param_r = 1.742; */

/* param_n = 14; */
/* param_r = sqrt(2);*/
/* param_r = 4.0 / 3.0; */
/* param_r = 5.0 / 4.0; */

/* param_n = 10; */
/* param_r = 1.545; */

zboxminx = -0.007;
zboxmaxx = 0.007;
zboxminy = -0.0064;
zboxmaxy = 0.0064;

min_thresh = 10^(-20);

wheight = sqrt(param_r^2 - 1);
wwidth = param_r - 1;

rm2d(t) = [cos(t), sin(t), 0; -1 * sin(t), cos(t), 0; 0, 0, 1];

tm2d(x, y) = [1, 0, x; 0, 1, y; 0, 0, 1];

/* twista = tm2d(-1, 0) * rm2d((2 * Pi) / param_n) * tm2d(1, 0); */

/* twistb = tm2d(1, 0) * rm2d((2 * Pi) / param_n) * tm2d(-1, 0); */

dist(p1, p2) = sqrt((p2[1] -  p1[1])^2 + (p2[2] - p1[2])^2);


rpoint() = [(random(1.0) * (2 * wwidth)) - wwidth, (random(1.0) * (2 * wheight)) - wheight, 1]~;

/* rpoint() = [(random(1.0) * 0.2) - 0.15, (random(1.0) * 0.15) -0.05, 1]~; */

inzbox(p) = if(p[1] >= zboxminx, if(p[1] <= zboxmaxx, if(p[2] >= zboxminy, if(p[2] <= zboxmaxy, 1, 0), 0), 0), 0);

incirca(param_r, p) = {if(dist(p, [-1, 0, 1]~) <= param_r, 1, 0);};
incircb(param_r, p) = {if(dist(p, [1, 0, 1]~) <= param_r, 1, 0);};

inwedge(param_r, p) = {if(incirca(param_r, p) == 1, if(incircb(param_r, p) == 1, 1, 0);, 0);}

goodrpoint() = {
found = 0;
while(found == 0, point = rpoint(); if(inwedge(point) == 1, found = 1, ));
point;
};

pointcycle(x, y) = {
count = 0;
startpoint = [x, y, 1]~;
point = startpoint;
printf("%.10f\n",point);
until(dist(point, startpoint) < min_thresh, if(incirca(point) == 1, point=twista*point; printf("%.10f\n",point); count = count + 1;,); if(incircb(point) == 1, point=twistb*point; printf("%.10f\n",point); count = count + 1;,););
printf("count = %d\n", count);
};

pointcycleorder(p) = {
count = 0;
startpoint = p;
point = startpoint;
until(dist(point, startpoint) < min_thresh, if(incirca(point) == 1, point=twista*point; count = count + 1;,); if(incircb(point) == 1, point=twistb*point; count = count + 1;,););
printf("%.10f %.10f %d\n", startpoint[1], startpoint[2] count);
};


pointcycleorderab(p, an, bn, l) = {
count = 0;
startpoint = p;
point = startpoint;
twistan = twista^an;
twistbn = twistb^bn;
until(dist(point, startpoint) < min_thresh, if(incirca(point) == 1, point=twistan*point; count = count + 1;,); if(incircb(point) == 1, point=twistbn*point; count = count + 1;,); if(count > l, break(),););
if(count > l, printf("# %.10f %.10f HIT LIMIT %d\n", startpoint[1], startpoint[2], l), printf("%.10f %.10f %d\n", startpoint[1], startpoint[2], count););
};

pointcycleorderabfast_order(p,an,bn)=count=0;pcount=0;startpoint=p;point=startpoint;twistan=twista^an;twistbn=twistb^bn;until(dist(point,startpoint)<min_thresh,if(incirca(point)==1,point=twistan*point;count=count+1;,);if(incircb(point)==1,point=twistbn*point;count=count+1;,);pcount=pcount+1;if(pcount>1000000,pcount=0;printf("Count so far: %d\n",count);,););printf("%.10f %.10f %d\n",startpoint[1],startpoint[2],count);


pointcycleorderabfast(param_n, param_r, p, an, bn, l) = {
count = 0;
startpoint = p;
point = startpoint;

twista = tm2d(-1, 0) * rm2d((2 * Pi) / param_n) * tm2d(1, 0);
twistb = tm2d(1, 0) * rm2d((2 * Pi) / param_n) * tm2d(-1, 0);

twistan = twista^an;
twistbn = twistb^bn;

wpoints = List();
until(dist(point, startpoint) < min_thresh,

if(incirca(param_r, point) == 1, point=twistan*point; count = count + 1; /*if(inwedge(point) == 1, listput(wpoints,point); ,);*/,);

if(incircb(param_r, point) == 1, point=twistbn*point; count = count + 1; if(inwedge(param_r, point) == 1, listput(wpoints,point);,);,);

if(count > l, break(),);

);

if(count > l, printf("# %.10f %.10f HIT LIMIT %d\n", startpoint[1], startpoint[2], l), wc = length(wpoints); for(x = 1, wc, printf("%.10f %.10f %d\n", wpoints[x][1], wpoints[x][2], count)););
};


pointcycleorderabfastzbox(p, an, bn, l) = {
count = 0;
startpoint = p;
point = startpoint;
twistan = twista^an;
twistbn = twistb^bn;
wpoints = List();
until(dist(point, startpoint) < min_thresh,

if(incirca(point) == 1, point=twistan*point; count = count + 1; /*if(inwedge(point) == 1, listput(wpoints,point); ,);*/,);

if(incircb(point) == 1, point=twistbn*point; count = count + 1; if(inzbox(point) == 1, if(inwedge(point) == 1, listput(wpoints,point);,);,);,);

if(count > l, break(),);

);

if(count > l, printf("# %.10f %.10f HIT LIMIT %d\n", startpoint[1], startpoint[2], l), wc = length(wpoints); for(x = 1, wc, printf("%.10f %.10f %d\n", wpoints[x][1], wpoints[x][2], count)););
};


pointcycle_ab_delta(param_n, param_r, p, an, bn, l) = {
counta = 0;
countb = 0;
startpoint = p;
point = startpoint;

twista = tm2d(-1, 0) * rm2d((2 * Pi) / param_n) * tm2d(1, 0);
twistb = tm2d(1, 0) * rm2d((2 * Pi) / param_n) * tm2d(-1, 0);

twistan = twista^an;
twistbn = twistb^bn;

wpoints = List();
until(dist(point, startpoint) < min_thresh,

if(incirca(param_r, point) == 1, point=twistan*point; counta = counta + 1; /*if(inwedge(point) == 1, listput(wpoints,point); ,);*/,);

if(incircb(param_r, point) == 1, point=twistbn*point; countb = countb + 1; if(inwedge(param_r, point) == 1, listput(wpoints,point);,);,);

if(counta + countb > l, break(),);

);

if(counta + countb > l, printf("# %.10f %.10f HIT LIMIT %d\n", startpoint[1], startpoint[2], l), wc = length(wpoints); for(x = 1, wc, printf("%.10f %.10f %d %d\n", wpoints[x][1], wpoints[x][2], counta, countb)););
};


/*for(x = 0, 10, pointcycleordera5b1(goodrpoint(), 5, 1));*/
/*for(x = 1, 100000, pointcycleorderabfast(goodrpoint(), 11, 1, 800));*/

xmin = -0.0106873;
xmax = 0.000116670;
ymin = 0.089184;
ymax = 0.11;
xres = 200;
yres = 200;
x = xmin;
/*while(x < xmax, y = ymin; while(y < ymax, pointcycleorderabfast([x, y, 1]~, 11, 1, 100000); y = y + ((ymax - ymin) / yres);); x = x + ((xmax - xmin) / xres)); */

/* \q */
