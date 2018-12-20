function [ CP ] = cpnew(protein, carbohydrate, fiber, fat, ash,temp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
cp1 = 2.0082 + (1.2089e-3)*temp - (1.3129e-6)*temp^2;
cp2 = 1.9842 + (1.4733e-3)*temp - (4.8008e-6)*temp^2;
cp3 = 1.5488 + (1.9625e-3)*temp - (5.9399e-6)*temp^2;
cp4 = 1.8459 + (1.8306e-3)*temp - (4.6509e-6)*temp^2;
cp5 = 1.0926 + (1.8896e-3)*temp - (3.6817e-6)*temp^2;

CP = protein*cp1 + carbohydrate*cp3 + fat*cp2 + fiber*cp4 + ash*cp5;
end

