#include<stdio.h>
#include<iostream>
#include<utility>
#include<cmath>
#include<opencv2/opencv.hpp>
#include<time.h>
using namespace std;
using namespace cv;
#define PI 3.141592653589793
#define length 2.5
#define velocity 90
#define time_del 0.1


typedef struct{
	bool val;
	pair<float,float> p;

} interStruct;

typedef struct{
	float delta;
	float orientation;
	pair<float,float> new_location;
} new_locStruct;

// bool,arr_contour,arr_orientation,delta_return
typedef struct{
	bool path_ret;
	vector<pair<float,float> > path;
	vector<float> orientation;
	vector<float> delta;

} pathStruct;

class purepursuit
{

public:
	pair<pair<float,float>,pair<float,float> > list_pts;
	float orientation;
	float delta_prev;
	
	purepursuit(pair<pair<float,float>,pair<float,float> > list_pts,float orientation,float delta_prev)
	{
		this->list_pts = list_pts;
		this->orientation = orientation;
		this->delta_prev = delta_prev;
	}


	float Distance(pair<float,float> pt1, pair<float,float> pt2){
		float value = sqrt((pt2.second - pt1.second)*(pt2.second - pt1.second) + (pt2.first - pt1.first)*(pt2.first - pt1.first));
		return value;
	}

	interStruct intersection_pt(pair<float,float> start,pair<float,float> end,pair<float,float> pt,float dist_look){

		interStruct ret;
		float lineDist = abs((end.second-start.second)*pt.first - (end.first-start.first)*pt.second + end.first*start.second - end.second*start.first)/sqrt((end.second-start.second)*(end.second-start.second)+(end.first-start.first)*(end.first-start.first));
		if(lineDist > dist_look){// pt_check(pt,make_pair(start,end)) == false){
			// cout<<"lineDist > dist_look"<<endl;
			ret.val = false;
			ret.p = make_pair(NAN,NAN);
			return ret;
		}

		else{
			ret.val = true;
			pair<float,float> I1;
			pair<float,float> I2;

			if(start.first == end.first){
				I1.first = start.first;
				I1.second = pt.second + sqrt(dist_look*dist_look - (start.first - pt.first)*(start.first - pt.first));

				I2.first = start.first;
				I2.second = pt.second - sqrt(dist_look*dist_look - (start.first - pt.first)*(start.first - pt.first));

			}
			else{
				float c = (start.second*end.first - start.first*end.second)/(end.first - start.first);
				float m = (end.second - start.second)/(end.first - start.first);
				float p = pt.first;
				float q = pt.second;
				float r = dist_look;

				float A = m*m + 1;			
				float B = 2*(m*c - m*q - p);
				float C = q*q - r*r + p*p - 2*c*q +c*c;
			
				I1.first = (-B + sqrt(B*B - 4*A*C))/(2*A);
				I1.second = m*(-B + sqrt(B*B - 4*A*C))/(2*A) + c;

				I2.first = (-B - sqrt(B*B - 4*A*C))/(2*A);
				I2.second = m*(-B - sqrt(B*B - 4*A*C))/(2*A) + c;
			}

			if(this->Distance(I1,end) > this->Distance(I2,end)){ret.p = I2;}
			else{ret.p = I1;}


			// if(pt_check(I1,make_pair(start,end)) == true && pt_check(I2,make_pair(start,end)) == true){
			// }
			// else if(pt_check(I1,make_pair(start,end)) == false && pt_check(I2,make_pair(start,end)) == false){
			// 	ret.val = false;
			// 	ret.p = make_pair(NAN,NAN);
			// }
			// else{
			// 	ret.p.first = pt_check(I1,make_pair(start,end))*I1.first + pt_check(I2,make_pair(start,end))*I2.first;
			// 	ret.p.second = pt_check(I1,make_pair(start,end))*I1.second + pt_check(I2,make_pair(start,end))*I2.second;
			// }

			return ret;
		}

	}

	float eta_fn(pair<float,float> pt,pair<float,float> pt_path,float orientation){

		float x_bot_vec = length*cos(orientation*PI/180.0);
		float y_bot_vec = length*sin(orientation*PI/180.0);
		float x_pt_vec = pt_path.first - pt.first;
		float y_pt_vec = pt_path.second - pt.second;

		float argument = ((x_bot_vec*x_pt_vec) + (y_bot_vec*y_pt_vec))/sqrt((x_bot_vec*x_bot_vec + y_bot_vec*y_bot_vec)*(x_pt_vec*x_pt_vec+y_pt_vec*y_pt_vec));
		if(abs(argument) > 1){argument = argument/abs(argument);}

		float eta = acos(argument)*180/PI;
		float k = y_bot_vec*x_pt_vec-x_bot_vec*y_pt_vec;

		if(k < 0){eta *= -1;}
		
		return eta;
	}

	float Delta_Deconv(float delta){
		if(delta >= 30){delta = delta - 30;}
		else{delta = 330 + delta;}

		return delta;
	}

	float Delta_conv(float delta){
		if(delta >= 330){delta = delta - 330;}
		else if(delta < 30){delta = 30 + delta;}

		return delta;
	}

	bool loopCheck(pair< pair<float,float>,pair<float,float> > list_pts,pair<float,float> prev_pt,pair<float,float> current_pt){

		float x1 = list_pts.first.first;
		float y1 = list_pts.first.second;
		float x2 = list_pts.second.first;
		float y2 = list_pts.second.first;

		float m1;
		if(list_pts.second.first != list_pts.first.first){
			m1 = (list_pts.second.second - list_pts.first.second)/(list_pts.second.first - list_pts.first.first);
		}
		else{m1 = NAN;}


		float m2;
		if(prev_pt.first != current_pt.first){
			m2 = (prev_pt.second - current_pt.second)/(prev_pt.first - current_pt.first);
		}
		else{m2 = NAN;}


		if(m1 != m2 && m1 != NAN && m2 != NAN){
			pair<float,float> pt;
			pt.first = ((y1-y2) - (m1*x1-m2*x2))/(m2-m1);
			pt.second = ((m2*y1-m1*y2) + m1*m2*(x2-x1))/(m2-m1);

			float PA = this->Distance(pt,list_pts.first);
			float AB = this->Distance(list_pts.first,list_pts.second);
		
			float PC;
			if(this->Distance(pt,prev_pt) > this->Distance(pt,current_pt)){
				PC = this->Distance(pt,prev_pt);
			}
			else{PC = this->Distance(pt,current_pt);}

			float CD = this->Distance(prev_pt,current_pt);

			if(PA/AB > 0.8 && PC/CD < 1){return true;}
			else{return false;}
		}

		else{return false;}
	}

	float mod(float val,float div){

		while(val < 0){
			val = val + div;
		}

		while(val > div){
			val = val - div;
		}
		return val;
	}

	new_locStruct new_loc(float eta,float dist_look_inter,pair<float,float> pt_old,float bot_orientation,float delta_prev)
	{
		new_locStruct next_loc;
		float delta_dot_max = 20;
		float Td = 0.05;
		float delta_next = 0;
		float theta_dot = 0;
		float time1;
		delta_prev = this->Delta_conv(delta_prev);

		float x_bot = pt_old.first;
		float y_bot = pt_old.second;

		float pt_new_x = 0;
		float pt_new_y = 0;
		float bot_orientation_new = 0;

		float delta = -(atan(2*length*sin(eta*PI/180)/dist_look_inter))*180/PI;
		delta = this->mod(delta,360);
		
		if(delta > 180 && delta < 330){delta = 330;}
		else if(delta < 180 && delta > 30){delta = 30;}
		delta = this->Delta_conv(delta);
			
		float delta_dot = (delta - delta_prev)/Td;

		if(delta_dot > delta_dot_max){delta_dot = delta_dot_max;}
		else if(delta_dot < -delta_dot_max){delta_dot = -delta_dot_max;}


		if(delta_dot == 0){time1 = 0;}
		else{time1 = (delta - delta_prev)/delta_dot;}

		if(time1 >= time_del){
			// cout<<"if"<<endl;
			for(int i = 0;i < 10;i++){
				delta_next = delta_prev + delta_dot*(i*time_del/10);
				
				if(delta_next > 60){delta_next = 60;}
				else if(delta_next < 0){delta_next = 0;}

				theta_dot = velocity*tan(this->Delta_Deconv(delta_next)*PI/180)/length;
				bot_orientation = this->mod(bot_orientation,360);

				if(theta_dot != 0){
					pt_new_x = x_bot + 2*velocity*cos(bot_orientation*PI/180 + theta_dot*time_del*PI/3600)*(sin(theta_dot*time_del*PI/3600))/theta_dot;
					pt_new_y = y_bot + 2*velocity*sin(bot_orientation*PI/180 + theta_dot*time_del*PI/3600)*(sin(theta_dot*time_del*PI/3600))/theta_dot;
					bot_orientation_new = bot_orientation + theta_dot*time_del/10;
				}
				else{
					pt_new_x = x_bot + velocity*cos(bot_orientation*PI/180)*time_del/10;
					pt_new_y = y_bot + velocity*sin(bot_orientation*PI/180)*time_del/10;
					bot_orientation_new = bot_orientation;
				}

				x_bot = pt_new_x;
				y_bot = pt_new_y;
				bot_orientation = bot_orientation_new;
				bot_orientation = this->mod(bot_orientation,360);

			}

		}

		else{
			for(int i = 0;i < 10 && time1 > 0;i++){
				delta_next = delta_prev + delta_dot*(i*time1/10);

				if(delta_next > 60){delta_next = 60;}
				else if(delta_next < 0){delta_next = 0;}

				theta_dot = velocity*tan(this->Delta_Deconv(delta_next)*PI/180)/length;
				bot_orientation = this->mod(bot_orientation,360);
				
				if(theta_dot != 0){	
					pt_new_x = x_bot + 2*velocity*cos(bot_orientation*PI/180 + theta_dot*time1*PI/3600)*(sin(theta_dot*time1*PI/3600))/theta_dot;
					pt_new_y = y_bot + 2*velocity*sin(bot_orientation*PI/180 + theta_dot*time1*PI/3600)*(sin(theta_dot*time1*PI/3600))/theta_dot;
					bot_orientation_new = bot_orientation + theta_dot*time1/10;
				}
				else{
					pt_new_x = x_bot + velocity*cos(bot_orientation*PI/180)*time1/10;
					pt_new_y = y_bot + velocity*sin(bot_orientation*PI/180)*time1/10;
					bot_orientation_new = bot_orientation;
				}

				x_bot = pt_new_x;
				y_bot = pt_new_y;
				bot_orientation = bot_orientation_new;
				bot_orientation = this->mod(bot_orientation,360);
			}


			theta_dot = velocity*tan(this->Delta_Deconv(delta)*PI/180)/length;
			bot_orientation = this->mod(bot_orientation,360);
			if(theta_dot != 0){
				pt_new_x = x_bot + 2*velocity*cos(bot_orientation*PI/180 + theta_dot*(time_del - time1)*PI/360)*(sin(theta_dot*(time_del - time1)*PI/360))/theta_dot;
				pt_new_y = y_bot + 2*velocity*sin(bot_orientation*PI/180 + theta_dot*(time_del - time1)*PI/360)*(sin(theta_dot*(time_del - time1)*PI/360))/theta_dot;
				bot_orientation_new = bot_orientation + theta_dot*(time_del - time1);
			}	
			else{
				pt_new_x = x_bot + velocity*cos(bot_orientation*PI/180)*(time_del - time1);			
				pt_new_y = y_bot + velocity*sin(bot_orientation*PI/180)*(time_del - time1);				
				bot_orientation_new = bot_orientation;
			}
			x_bot = pt_new_x;
			y_bot = pt_new_y;
			bot_orientation = bot_orientation_new;
			bot_orientation = this->mod(bot_orientation,360);
		}

		delta_next = this->Delta_Deconv(delta_next);

		next_loc.delta = delta_next;
		next_loc.orientation = bot_orientation;
		next_loc.new_location.first = x_bot;
		next_loc.new_location.second = y_bot; 
		return next_loc;
	}

	pathStruct path()
	{
		pair<pair<float,float>,pair<float,float> > list_pts = this->list_pts;
		float orientation = this->orientation;
		float delta_prev = this->delta_prev;


		pair<float,float> pt_bot;
		float x_bot = list_pts.first.first;
		float y_bot = list_pts.first.second;
		float bot_orientation = orientation;
		float bot_orientation_new = 0;
		float prev_index = 0;
		float eta = 0;
		float dist_look = 15;
		float delta;
		
		vector<pair<float,float> >  arr_contour;
		vector<float> 				arr_orientation;
		vector<float> 				delta_return;

		arr_contour.push_back(list_pts.first);
		arr_orientation.push_back(bot_orientation);
		pair<float,float> pt_path;
		pair<float,float> pt_new;

		float pt_new_x = 0;
		float pt_new_y = 0;
		bool val_sample = 0;
		bool val_loopcheck = 0;
		bool val_return = false;
		bool val_lookAhead = false;

		while(true)
		{
			pt_bot = make_pair(x_bot,y_bot);
			if(val_lookAhead == false){
				interStruct intersect = this->intersection_pt(list_pts.first,list_pts.second,pt_bot,dist_look);
				val_sample = intersect.val;
				pt_path = intersect.p;

				if(val_sample == false){
					pathStruct ret;
					ret.path_ret = false;
					ret.path.push_back(make_pair(NAN,NAN));
					ret.orientation.push_back(NAN);
					ret.delta.push_back(NAN);
					return ret;
				}
			}

			else{
				pt_path = list_pts.second;
				dist_look = this->Distance(pt_path,pt_bot)/2;
			}	
			eta = this->eta_fn(pt_bot,pt_path,bot_orientation);

			// cout<<"pt_bot : ["<<pt_bot.first<<','<<pt_bot.second<<"] ,pt_path : ["<<pt_path.first<<','<<pt_path.second<<"],orien : "<<bot_orientation<<endl;
			// cout<<"Delta : "<<delta<<", eta : "<<eta<<endl;
			// cout<<endl;
			new_locStruct new_loc_details = this->new_loc(eta,dist_look,pt_bot,bot_orientation,delta_prev);
			delta = new_loc_details.delta;
			pt_new = new_loc_details.new_location;
			bot_orientation_new = new_loc_details.orientation;

			if(this->Distance(pt_bot,list_pts.second) < 10){
				val_lookAhead = true;
				val_loopcheck = this->loopCheck(list_pts,pt_bot,pt_new);
				if(val_loopcheck == true){val_return = true;}
			}

			pt_new_x = pt_new.first;
			pt_new_y = pt_new.second;
			
			delta_return.push_back(delta);
			arr_contour.push_back(pt_new);
			arr_orientation.push_back(bot_orientation_new);

			if(val_return == true || this->Distance(pt_new,list_pts.second) < 5){
				pathStruct ret;
				ret.path_ret = true;
				ret.path = arr_contour;
				ret.orientation = arr_orientation;
				ret.delta = delta_return;

				return ret;
			}
			
			x_bot = pt_new_x;
			y_bot = pt_new_y;
			pt_bot = pt_new;
			bot_orientation = bot_orientation_new;
			delta_prev = delta;

		}
		pathStruct some;
		return some;
	}

};

void MyLine( Mat img, Point start, Point end,unsigned int val)
{
  int thickness = 1;
  int lineType = 8;
  line( img,start,end,Scalar(val,val,val),thickness,lineType );
}

int main(){

	float orientation = 90;
	float delta = 0;
	float dist_look = 15;
	
	pair< pair<float,float>,pair<float,float> > list_pts = make_pair(make_pair(100,100),make_pair(295,170));	
	purepursuit p(list_pts,orientation,delta);

	clock_t t1,t2;
	t1 = clock();
	pathStruct path_return = p.path();
	t2 = clock();
	cout<<(float)(t2-t1)/CLOCKS_PER_SEC<<endl;

	cout<<path_return.path.size()<<endl;
	vector<pair<float,float> > arr_contour = path_return.path;
	Mat image = Mat::zeros( 500,500,CV_8UC3);
	MyLine(image,Point(list_pts.first.first,list_pts.first.second),Point(list_pts.second.first,list_pts.second.second),125);
	namedWindow("image");
	for(int i = 0;i<path_return.path.size()-1;i++){
		MyLine(image,Point(abs(arr_contour[i].first),abs(arr_contour[i].second)),Point(abs(arr_contour[i+1].first),abs(arr_contour[i+1].second)),255);
		imshow("image",image);
		waitKey(1);		
	}
	waitKey(0);

	return 0;

}

