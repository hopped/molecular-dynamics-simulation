/***
Copyright 2013 University of Stuttgart  
Author: Anthony Sulistio (HLRS)
A header file for converting pov to png images
***/  

#declare White = rgb <1,1,1>;
#declare LightGray = White*0.8;

// taken from http://www.f-lohmueller.de/pov_tut/basic/povtuto4.htm#camera
/****
camera {  //ultra_wide_angle
         orthographic angle 10
         location<0,0,-10> // standpoint of the viewer
         look_at <2,2,-5> // where to look at
         //right x*image_width/image_height  // aspect ratio
         angle 70 // camera angle
       } 
******/

/*******/ // for 100k molecules
camera{ orthographic angle 110
        location <3,3,-3>*10
        look_at <5,5,-5> // where to look at
        right x*image_width/image_height
        translate <35,35,0>
      }
/*****/

/*****
// for 500k and 600k molecules
camera{ orthographic angle 110
        location <3,3,-3>*10
        look_at <-8,-8,8> // where to look at
        right x*image_width/image_height
        translate <88,88,-10>
      }
******/

//background { color rgb <1, 1, 1> }   // white
background { color LightGray }

//light_source { <1,1,-3> color 1 }
//light_source { <-10,-20,-20> color White }
light_source{ <3,3,0> color rgb <1,1,1>
              parallel
              point_at<1, 0, 0> 
            }

