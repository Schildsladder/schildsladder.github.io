function draw() {
    const canvas = document.getElementById("graph");
    if (canvas.getContext) {
      const ctx = canvas.getContext("2d");

      ctx.fillStyle = "rgb(200, 0, 0)";
      ctx.fillRect(10, 10, 50, 50);

      ctx.fillStyle = "rgba(0, 0, 200, 0.5)";
      ctx.fillRect(30, 30, 50, 50);
    }
  }
draw();



var canvas = document.getElementById("graph");
var context = canvas.getContext('2d');

// source for this code: https://riptutorial.com/html5-canvas/example/11659/detecting-mouse-position-on-the-canvas

canvas.addEventListener("mousemove", function(e) { 
    var cRect = canvas.getBoundingClientRect();
    var canvasX = Math.round(e.clientX - cRect.left);
    var canvasY = Math.round(e.clientY - cRect.top);
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.fillText("X: "+canvasX+", Y: "+canvasY, 10, 20);
});





ctx.beginPath();
ctx.moveTo(75, 50);
ctx.lineTo(100, 75);
ctx.lineTo(100, 25);
ctx.fill();

ctx.beginPath();
ctx.arc(75, 75, 50, 0, Math.PI * 2, true); // Outer circle
ctx.moveTo(110, 75);
ctx.arc(75, 75, 35, 0, Math.PI, false); // Mouth (clockwise)
ctx.stroke();

ctx.lineTo(45, 125);
ctx.closePath();
ctx.stroke();


function draw() {
    const ctx = document.getElementById("canvas").getContext("2d");
    for (let i = 0; i < 6; i++) {
      for (let j = 0; j < 6; j++) {
        ctx.fillStyle = `rgb(${Math.floor(255 - 42.5 * i)}, ${Math.floor(
          255 - 42.5 * j
        )}, 0)`;
        ctx.fillRect(j * 25, i * 25, 25, 25);
      }
    }
  }
  

Instance methods
    arc()
    arcTo()
    beginPath()
    bezierCurveTo()
    clearRect()
    clip()
    closePath()
    createConicGradient()
    createImageData()
    createLinearGradient()
    createPattern()
    createRadialGradient()
    drawFocusIfNeeded()
    drawImage()
    ellipse()
    fill()
    fillRect()
    fillText()
    getContextAttributes()
    getImageData()
    getLineDash()
    getTransform()
    Experimental
    isContextLost()
    isPointInPath()
    isPointInStroke()
    lineTo()
    measureText()
    moveTo()
    putImageData()
    quadraticCurveTo()
    rect()
    Experimental
    reset()
    resetTransform()
    restore()
    rotate()
    roundRect()
    save()
    scale()
    Experimental
    scrollPathIntoView()
    setLineDash()
    setTransform()
    stroke()
    strokeRect()
    strokeText()
    transform()
    translate()

    Instance properties

    canvas
    direction

    fillStyle
    ctx.fillStyle = "orange";
    ctx.fillStyle = "#FFA500";
    ctx.fillStyle = "rgb(255, 165, 0)";
    ctx.fillStyle = "rgba(255, 165, 0, 1)";

    filter
    font
    fontKerning
    Experimental
    fontStretch
    Experimental
    fontVariantCaps
    globalAlpha  global transparancy ctx.globalAlpha = 0.2;
    globalCompositeOperation
    imageSmoothingEnabled
    imageSmoothingQuality
    Experimental
    letterSpacing
    lineCap
    lineDashOffset
    lineJoin
    lineWidth ex =>  ctx.lineWidth = 15;
    miterLimit
    shadowBlur
    shadowColor
    shadowOffsetX
    shadowOffsetY

    strokeStyle
    ctx.strokeStyle = "rgba(255, 0, 0, 0.5)";  

    textAlign
    textBaseline
    Experimental
    textRendering
    Experimental
    wordSpacing