#// Apperture diffraction is the Fourier Transform of the shape of the apperture.
#// Turns out you can derive some close forms for simple apperture shapes.
#// With mathematica I got 3,4,6 bladed (https://physics.stackexchange.com/questions/9899/how-does-fraunhofer-diffraction-depend-on-the-orientation-of-the-sides-of-a-lens)
#// Circle apperture is deffined by Airy Diffraction.
#// 
#// Be sure to do the math in full spectral space and integrate down to RGB for better results
#
#// supported apperture blades 3-4-6-Circle
#const int NumberOfAppertureBlades = 4;
#
#
#// Spectrum to xyz approx function from "Simple Analytic Approximations to the CIE XYZ Color Matching Functions"
#// http://jcgt.org/published/0002/02/01/paper.pdf
#//Inputs:  Wavelength in nanometers
#float xFit_1931(float wave)
#{
#    float t1 = (wave - 442.0f)*((wave < 442.0f) ? 0.0624f : 0.0374f);
#    float t2 = (wave - 599.8f)*((wave < 599.8f) ? 0.0264f : 0.0323f);
#    float t3 = (wave - 501.1f)*((wave < 501.1f) ? 0.0490f : 0.0382f);
#    return 0.362f*exp(-0.5f*t1*t1) + 1.056f*exp(-0.5f*t2*t2) - 0.065f*exp(-0.5f*t3*t3);
#}
#float yFit_1931(float wave)
#{
#    float t1 = (wave - 568.8f)*((wave < 568.8f) ? 0.0213f : 0.0247f);
#    float t2 = (wave - 530.9f)*((wave < 530.9f) ? 0.0613f : 0.0322f);
#    return 0.821f*exp(-0.5f*t1*t1) + 0.286f*exp(-0.5f*t2*t2);
#}
#float zFit_1931(float wave)
#{
#    float t1 = (wave - 437.0f)*((wave < 437.0f) ? 0.0845f : 0.0278f);
#    float t2 = (wave - 459.0f)*((wave < 459.0f) ? 0.0385f : 0.0725f);
#
#    return 1.217f*exp(-0.5f*t1*t1) + 0.681f*exp(-0.5f*t2*t2);
#}
#
#// RGB-XYZ Matrix Calculator
#// http://www.russellcottrell.com/photo/matrixCalculator.htm
#// Based on equations found here:
#// http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
#// And Rec. 2020 values found here:
#// https://en.wikipedia.org/wiki/Rec._2020
#// https://en.wikipedia.org/wiki/Rec._709
#// https://en.wikipedia.org/wiki/SRGB
#vec3 XYZtosRGB(vec3 XYZ)
#{
#    vec3 rgb;
#    rgb.x = XYZ.x *  3.2409699f + XYZ.y * -1.5373832f + XYZ.z * -0.4986108f;
#    rgb.y = XYZ.x * -0.9692436f + XYZ.y *  1.8759675f + XYZ.z *  0.0415551f;
#    rgb.z = XYZ.x *  0.0556301f + XYZ.y * -0.2039770f + XYZ.z *  1.0569715f;
#    
#    return rgb;
#}
#
#// usefull functions
##define M_2PI 6.28318530718
#float sinc(float v)
#{
#    float res = 1.0;
#    if (abs(v) > 0.0001)
#        res = sin(v) / v;
#
#    return res;
#}
#
##define SPECTRAL_SAMPLES 32
#
#
#
#void mainImage( out vec4 fragColor, in vec2 fragCoord )
#{
#    // Normalized pixel coordinates (from 0 to 1)
#    vec2 dd = ( fragCoord - .5*iResolution.xy ) / iResolution.y;
#    float d = length(dd);
#
#    
#    // https://en.wikipedia.org/wiki/Airy_disk
#    // cone spacing in the human eye are around 2.5um, resulting of an apperture of 5.0um
#    const float Apperture = 2.0*2.5;
#    
#    float lamdaStart = 380.0f;
#    float lamdaEnd = 780.0f;
#
#    float dw = (lamdaEnd - lamdaStart) / float(SPECTRAL_SAMPLES);
#
#    // https://en.wikipedia.org/wiki/CIE_1931_color_space
#    // integrating Illuminant
#    vec3 XYZIrradiance = vec3(0.0);
#
#    for (int i = 0; i < SPECTRAL_SAMPLES; i++)
#    {
#        float w = lamdaStart + (float(i) + 0.5f)*dw;
#        
#        float k = M_2PI / w; // k = 2PI/lamda in nm-1
#        float x = d*k*Apperture*1000.0; // um * nm-1 = 10^3
#
#        float dx = dd.x*k*Apperture*1000.0;
#        float dy = dd.y*k*Apperture*1000.0;
#
#        float illuminant = 0.0;
#        
#        x *= 2.0;
#        dx *= 2.0;
#        dy *= 2.0;
#if (NumberOfAppertureBlades == 3)
#{
#    	// 3 Apperture Blades
#        float u = dx;
#        float u2 = u*u;
#        float u4 = u2*u2;
#
#        float v = dy;
#        float v2 = v*v;
#        illuminant = (3.0*u2 +
#            3.0*v2 + (u2 - 3.0*v2)*cos(3.0*u) -
#            2.0*u*(u + sqrt(3.0)*v)*cos(3.0 / 2.0*(u - sqrt(3.0)*v)) -
#            2.0*u2*cos(3.0 / 2.0*(u + sqrt(3.0)*v)) +
#            2.0*sqrt(3.0)*u*v*cos(3.0 / 2.0*(u + sqrt(3.0)*v)));
#        illuminant *= 3.0 / (2.0*3.14*3.14*u2*(u2 - 3.0*v2)*(u2 - 3.0*v2));
#        if ((length(u - sqrt(3.0)*v) < 0.001) || (length(-u - sqrt(3.0)*v) < 0.001))
#        {
#
#            float s = sin(3.0*u / 2.0);
#            illuminant = 3.0*s*s*s*s / (4.0*3.14*3.14*u4);
#        }
#}
#else if (NumberOfAppertureBlades == 4)
#{
#        // 4 Apperture Blades
#        illuminant = sinc(dx)*sinc(dx)*sinc(dy)*sinc(dy);
#}
#else if(NumberOfAppertureBlades == 6)
#{
#        // 6 Apperture Blades
#		illuminant = (2.0*sqrt(3.0)*dx*(cos(dx / 2.0)*cos(sqrt(3.0)*dy / 2.0) - cos(dx)) -
#                    6.0*dy*sin(dx / 2.0)*sin(sqrt(3.0)*dy / 2.0)) / (3.14*dx*(dx*dx - 3.0*dy*dy));
#        illuminant = illuminant*illuminant;
#}
#else
#{        
#        // Circle apperture
#        // (2BesselJ(1,x)/x)^2 ~ cos^2(x-2Pi/4)/x^3
#        illuminant = cos(x-2.0*3.14/4.0)* cos(x-2.0*3.14/4.0)/(x*x*x);
#    
#        // [Muhammad Taher Abuelma’atti] Trigonometric Approximations for some Bessel Functions
#   //     illuminant = 1.0/6.0*sin(x/2.0) + 1.0/6.0*sin(x) + sqrt(3.0)/6.0*sin(sqrt(3.0)*x/2.0);
#   //     illuminant = (2.0*illuminant/x)*(2.0*illuminant/x);
#            
#}
#        XYZIrradiance.x += illuminant * xFit_1931(w) * dw;
#        XYZIrradiance.y += illuminant * yFit_1931(w) * dw;
#        XYZIrradiance.z += illuminant * zFit_1931(w) * dw;
#    }
#
#    
#    // XYZ to sRGB
#    vec3 RGBIrradiance = XYZtosRGB(XYZIrradiance);
# 	// sRGB
#    float a = 0.05;
#    RGBIrradiance.r = RGBIrradiance.r <= 0.0031308 ? 12.92*RGBIrradiance.r : (1.0+a)*pow(RGBIrradiance.r,1.0/2.4)-a;
#    RGBIrradiance.g = RGBIrradiance.g <= 0.0031308 ? 12.92*RGBIrradiance.g : (1.0+a)*pow(RGBIrradiance.g,1.0/2.4)-a;
#    RGBIrradiance.b = RGBIrradiance.b <= 0.0031308 ? 12.92*RGBIrradiance.b : (1.0+a)*pow(RGBIrradiance.b,1.0/2.4)-a;
#    
#    
#     fragColor = vec4(RGBIrradiance,1.0);
#
#}

#python implementation of the above shader
import math

def sinc(x):
    if x == 0:
        return 1
    else:
        return math.sin(x)/x

NumberOfAppertureBlades = 4

def xFit_1931(w):
    if w < 440:
        return 0.0014
    elif w < 490:
        return (w-440.0)/(490.0-440.0)*0.0624 + (490.0-w)/(490.0-440.0)*0.0014
    elif w < 510:
        return 0.0624
    elif w < 580:
        return (w-510.0)/(580.0-510.0)*0.3724 + (580.0-w)/(580.0-510.0)*0.0624
    elif w < 645:
        return (w-580.0)/(645.0-580.0)*0.2124 + (645.0-w)/(645.0-580.0)*0.3724
    elif w < 781:
        return (w-645.0)/(781.0-645.0)*0.0524 + (781.0-w)/(781.0-645.0)*0.2124
    else:
        return 0.0014

def yFit_1931(w):
    if w < 440:
        return 0.0003
    elif w < 490:
        return (w-440.0)/(490.0-440.0)*0.0134 + (490.0-w)/(490.0-440.0)*0.0003
    elif w < 510:
        return 0.0134
    elif w < 580:
        return (w-510.0)/(580.0-510.0)*0.1074 + (580.0-w)/(580.0-510.0)*0.0134
    elif w < 645:
        return (w-580.0)/(645.0-580.0)*0.0744 + (645.0-w)/(645.0-580.0)*0.1074
    elif w < 781:
        return (w-645.0)/(781.0-645.0)*0.0144 + (781.0-w)/(781.0-645.0)*0.0744
    else:
        return 0.0003

def zFit_1931(w):
    if w < 440:
        return 0.0065
    elif w < 490:
        return (w-440.0)/(490.0-440.0)*0.1159 + (490.0-w)/(490.0-440.0)*0.0065
    elif w < 510:
        return 0.1159
    elif w < 580:
        return (w-510.0)/(580.0-510.0)*0.3736 + (580.0-w)/(580.0-510.0)*0.1159
    elif w < 645:
        return (w-580.0)/(645.0-580.0)*0.3362 + (645.0-w)/(645.0-580.0)*0.3736
    elif w < 781:
        return (w-645.0)/(781.0-645.0)*0.0776 + (781.0-w)/(781.0-645.0)*0.3362
    else:
        return 0.0065

def XYZtosRGB(XYZ):
    return [3.2404542*XYZ[0] - 1.5371385*XYZ[1] - 0.4985314*XYZ[2],
            -0.9692660*XYZ[0] + 1.8760108*XYZ[1] + 0.0415560*XYZ[2],
            0.0556434*XYZ[0] - 0.2040259*XYZ[1] + 1.0572252*XYZ[2]]

SPECTRAL_SAMPLES = 100

def illuminant(w):
    x = w/1000.0
    return math.cos(x-2.0*math.pi/4.0)* math.cos(x-2.0*math.pi/4.0)/(x*x*x)

def main(fragCoord):
    XYZIrradiance = [0.0, 0.0, 0.0]
    for i in range(SPECTRAL_SAMPLES):
        w = 380.0 + i*(780.0-380.0)/SPECTRAL_SAMPLES
        dw = (780.0-380.0)/SPECTRAL_SAMPLES
        illuminant_ = illuminant(w)
        XYZIrradiance[0] += illuminant_ * xFit_1931(w) * dw
        XYZIrradiance[1] += illuminant_ * yFit_1931(w) * dw
        XYZIrradiance[2] += illuminant_ * zFit_1931(w) * dw
    RGBIrradiance = XYZtosRGB(XYZIrradiance)
    a = 0.05
    RGBIrradiance[0] = RGBIrradiance[0] <= 0.0031308 and 12.92*RGBIrradiance[0] or (1.0+a)*pow(RGBIrradiance[0],1.0/2.4)-a
    RGBIrradiance[1] = RGBIrradiance[1] <= 0.0031308 and 12.92*RGBIrradiance[1] or (1.0+a)*pow(RGBIrradiance[1],1.0/2.4)-a
    RGBIrradiance[2] = RGBIrradiance[2] <= 0.0031308 and 12.92*RGBIrradiance[2] or (1.0+a)*pow(RGBIrradiance[2],1.0/2.4)-a
    return RGBIrradiance
