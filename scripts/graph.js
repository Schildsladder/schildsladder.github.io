const earth = new Image();
var a=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

var f=[  0.001492,  -0.001633,  0.001756,  -0.001764,  0.001406,  -0.002146,  0.001451,  -0.001758,
    0.001860,  -0.001931,  0.001239,  -0.002382,  0.001371,  -0.002041,  0.002075,  -0.002310,
    0.000854,  -0.002890,  0.001178,  -0.002777,  0.002602,  -0.003342,  -0.000245,  -0.004292,
    0.000546,  -0.005222,  0.004272,  -0.007229,  -0.004803,  -0.010506,  -0.003109,  -0.018372,
    0.012604,  -0.036407,  -0.050415,  -0.102356,  -0.112793,  -0.254639,  -0.749023,  -0.335693,
    -0.000000,  -7.050781,  -0.782715,  0.908203,  -0.265381,  0.133057,  0.053192,  0.020493,
    0.003008,  0.022995,  -0.010796,  0.009697,  0.007889,  0.003069,  0.002325,  0.002974,
    -0.001130,  0.001263,  0.003389,  0.000319,  0.001945,  -0.000193,  0.000603,  -0.000525,
    0.002312,  -0.000525,  0.001796,  -0.001079,  0.001130,  -0.001109,  0.001936,  -0.000885,
    0.001725,  -0.001404,  0.001342,  -0.001351,  0.001775,  -0.001073,  0.001687,  -0.001544]

var n=[-40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23,
    -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10,  -9,  -8,  -7,  -6,  -5,
     -4,  -3,  -2,  -1,   0,   1,   2,   3,   4,   5,   6,  7,   8,   9,  10,  11,  12,  13,
     14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39]
   
let t=0;
const fps=60;
var heartx=[]
var hearty=[]




function init(){
    earth.src = "images/site_icon.ico";
    window.requestAnimationFrame(draw);
}

function draw() {
    const canvas = document.getElementById("canvas");
    const ctx = canvas.getContext("2d");

    ctx.globalCompositeOperation = "destination-over";
    ctx.clearRect(0, 0, 400, 400); // clear canvas
    // ctx.fillStyle = "rgba(0, 0, 0, 0.4)";
    ctx.strokeStyle = "rgba(255, 4, 2, 0.8)";
    ctx.save();


    var theta=t/500;
    t++;
    if (t==500){
        t=0;
    }


    ctx.translate(200, 200);
    ctx.beginPath();
    ctx.moveTo(0,0);

    var imag=0;
    var real=0;
    for (j=0; j<80  ;j++){


        imag=imag+f[j]*Math.sin(2*Math.PI*n[j]*theta);
        real=real-f[j]*Math.cos(2*Math.PI*n[j]*theta);

        ctx.lineTo(imag*20, real*20)
    }

    ctx.stroke();   
    ctx.restore();
    ctx.save()

    heartx.push(imag)
    heartx=heartx.slice(-100,)

    hearty.push(real)
    hearty=hearty.slice(-100,)

    ctx.translate(200, 200);
    ctx.beginPath();
    ctx.moveTo(heartx[0]*20,hearty[0]*20);
    
    for (j=0;j<100;j++){
        ctx.lineTo(heartx[j]*20,hearty[j]*20)
    }
    ctx.stroke();   
    ctx.restore();
    ctx.save()

    setTimeout(()=> {window.requestAnimationFrame(draw);}, 1000/fps);
  }

init();
