//	animate.js
//
//	Core routines for animated JavaScript applet that starts when
//	page is loaded.

//	To specialise, override customAnimate.initApp1()
//	and customAnimate.draw1(t) in a script included after this.

customAnimate = {};
customAnimate.initApp1 = function () {};
customAnimate.draw1 = function (t) {};

window.addEventListener('load',initApp0,false);
window.addEventListener('pageshow',initApp0,false);
window.addEventListener('beforeunload',stopApp,false);
window.addEventListener('unload',stopApp,false);
wkAnim = ('webkitRequestAnimationFrame' in window);
msAnim = ('msRequestAnimationFrame' in window);
mozAnim = ('mozRequestAnimationFrame' in window);
npAnim = ('requestAnimationFrame' in window);
wkcAnim = ('webkitCancelAnimationFrame' in window);
mscAnim = ('msCancelAnimationFrame' in window);
mozcAnim = ('mozCancelAnimationFrame' in window);
npcAnim = ('cancelAnimationFrame' in window);

canvas=null;
ctx=null;
timeBase=-1;
previousTime=-1;
animationTime = 60;
requestID=null;
animRequests=-1;

animControls = null;
animPauseButton = null;
animRedoButton = null;
animHideReq = null;
appPaused = false;
animationPaused = false;

function initApp0()
{
if (animRequests < 0)
	{
	animRequests = 0;
	initApp();
	};
}

function initApp()
{
timeBase=-1;
if (canvas=document.getElementById('animate'))
	{
	customAnimate.initApp1();
	queueAnim();
	appPaused = false;
	animationPaused = false;
	};
}

function queueAnim()
{
if (animRequests==0)
	{
	if (npAnim) requestID=window.requestAnimationFrame(draw);
	else if (wkAnim) requestID=window.webkitRequestAnimationFrame(draw);
	else if (msAnim) requestID=window.msRequestAnimationFrame(draw);
	else if (mozAnim) requestID=window.mozRequestAnimationFrame(draw);
	else requestID=window.setTimeout(draw,animationTime);
	animRequests++;
	};
}

function stopApp()
{
if (requestID!=null)
	{
	if (npcAnim) window.cancelAnimationFrame(requestID);
	else if (wkcAnim) window.webkitCancelAnimationFrame(requestID);
	else if (mscAnim) window.msCancelAnimationFrame(requestID);
	else if (mozcAnim) window.mozCancelAnimationFrame(requestID);
	else window.clearTimeout(requestID);
	requestID=null;
	};
animRequests = -1;

if (animHideReq) window.clearTimeout(animHideReq);
hideControls();
}

function draw(tt)
{
animRequests--;
requestID=null;

var currentTime=Date.now();
if (timeBase<0)
	{
	timeBase=currentTime;
	}
else if (currentTime - previousTime < animationTime)
	{
	queueAnim();
	return;
	};

previousTime = currentTime;

customAnimate.draw1((currentTime-timeBase)/1000.0);

queueAnim();
}

function animationControls()
{
if (!canvas) return;

var canvasRect = canvas.getBoundingClientRect(), st;
if (!animControls)
	{
	animControls = document.createElement("fieldset");
	document.body.appendChild(animControls);
	st=animControls.style;
	st.position="absolute";
	st.zIndex=1;
	
	animPauseButton = document.createElement("button");
	animControls.appendChild(animPauseButton);
	animPauseButton.onclick = pauseClick;

	animRedoButton = document.createElement("button");
	animControls.appendChild(animRedoButton);
	animRedoButton.value = animRedoButton.innerHTML = "Redraw";
	animRedoButton.onclick = redoClick;
	};
	
st=animControls.style;
var x = (window.pageXOffset !== undefined)
	? window.pageXOffset
	: (document.documentElement || document.body.parentNode || document.body).scrollLeft;
var y = (window.pageYOffset !== undefined)
	? window.pageYOffset
	: (document.documentElement || document.body.parentNode || document.body).scrollTop;
st.top=Math.round(canvasRect.top+y)+"px";
st.left=Math.round(canvasRect.left+x)+"px";
	
labelPauseButton();
	
if (animControls) animControls.style.display="block";
if (animHideReq) window.clearTimeout(animHideReq);
animHideReq = window.setTimeout(hideControls,5000);
}

function hideControls()
{
if (animControls) animControls.style.display="none";
animHideReq = null;
}

function labelPauseButton()
{
if (animPauseButton) animPauseButton.value = animPauseButton.innerHTML = appPaused ? "Resume" : "Pause";
}

function pauseClick()
{
appPaused = !appPaused;
labelPauseButton();

if (appPaused)
	{
	animationControls();
	stopApp();
	}
else
	{
	animRequests = 0;
	queueAnim();
	if (animHideReq) window.clearTimeout(animHideReq);
	hideControls();
	};
}

function redoClick()
{
if (animHideReq) window.clearTimeout(animHideReq);
hideControls();
animRequests = 0;
initApp();
}

function pauseAnimation()
{
if (!animationPaused)
	{
	animationPaused=true;
	stopApp();
	};
}

function resumeAnimation()
{
if (animationPaused)
	{
	animationPaused=false;
	animRequests=0;
	queueAnim();
	};
}
