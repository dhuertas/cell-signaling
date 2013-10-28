// Simulation environment -------------------------------------------------- //
var settings;

var spaceSize;

var particles = [];

var xhrStream;

var spheresMap = {};

var sphereMaterial;

// Camera coordinates
var r, theta, phi;

var x, y, z;

var max;

var scene, camera;
var renderer;

var pointLight;

var xAxis, yAxis, zAxis;
var bottomSide, topSide, leftSide, rightSide;

// Variables to parse the received data so far
var start, end;

start = 0;
end = 0;

var view = {
	width: 0,
	height: 0
};

function initEnvironment() {

	view.width = document.getElementById("view").clientWidth;
	view.height = document.getElementById("view").clientHeight;

	max = Math.max(spaceSize.x, spaceSize.y, spaceSize.z);

	scene = new THREE.Scene();
	camera = new THREE.PerspectiveCamera( 45, view.width / view.height, 0.1, 5*max);

	r = 2*max;
	theta = 0;
	phi = 0;

	x = spaceSize.x/2 + r*Math.sin(theta)*Math.cos(phi);
	y = spaceSize.y/2 + r*Math.sin(theta)*Math.sin(phi);
	z = spaceSize.z/2 + r*Math.cos(theta);

	camera.position.x = x;
	camera.position.y = y;
	camera.position.z = z;

	renderer = new THREE.WebGLRenderer();

	renderer.setSize(view.width, view.height
		/* window.innerWidth, window.innerHeight */
	);

	document.getElementById("view").appendChild(renderer.domElement);

	var material = new THREE.LineBasicMaterial({ 
		color: 0x000000
	});

	// axis
	var geo = new THREE.Geometry();

	geo.vertices.push(new THREE.Vector3(-max*0.1, -max*0.1, -max*0.1));		
	geo.vertices.push(new THREE.Vector3(0, -max*0.1, -max*0.1));

	xAxis = new THREE.Line(geo, new THREE.LineBasicMaterial({ color: 0xff0000 }));

	var geo = new THREE.Geometry();

	geo.vertices.push(new THREE.Vector3(-max*0.1, -max*0.1, -max*0.1));		
	geo.vertices.push(new THREE.Vector3(-max*0.1, 0, -max*0.1));

	yAxis = new THREE.Line(geo, new THREE.LineBasicMaterial({ color: 0x00ff00 }));

	var geo = new THREE.Geometry();

	geo.vertices.push(new THREE.Vector3(-max*0.1, -max*0.1, -max*0.1));		
	geo.vertices.push(new THREE.Vector3(-max*0.1, -max*0.1, 0));

	zAxis = new THREE.Line(geo, new THREE.LineBasicMaterial({ color: 0x0000ff }));

	scene.add(xAxis);
	scene.add(yAxis);
	scene.add(zAxis);

	// sim space bottom side
	var geometry = new THREE.Geometry();

	geometry.vertices.push(new THREE.Vector3(          0,           0, 0));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x,           0, 0));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x, spaceSize.y, 0));
	geometry.vertices.push(new THREE.Vector3(          0, spaceSize.y, 0));
	geometry.vertices.push(new THREE.Vector3(          0,           0, 0));

	bottomSide = new THREE.Line(geometry, material);

	scene.add(bottomSide);

	// sim space top side
	var geometry = new THREE.Geometry();

	geometry.vertices.push(new THREE.Vector3(          0,           0, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x,           0, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x, spaceSize.y, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(          0, spaceSize.y, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(          0,           0, spaceSize.z));

	topSide = new THREE.Line(geometry, material);

	scene.add(topSide);

	// sim space left side
	var geometry = new THREE.Geometry();

	geometry.vertices.push(new THREE.Vector3(          0, 0,           0));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x, 0,           0));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x, 0, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(          0, 0, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(          0, 0,           0));

	leftSide = new THREE.Line(geometry, material);

	scene.add(leftSide);

	// sim space right side
	var geometry = new THREE.Geometry();

	geometry.vertices.push(new THREE.Vector3(          0, spaceSize.y,           0));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x, spaceSize.y,           0));
	geometry.vertices.push(new THREE.Vector3(spaceSize.x, spaceSize.y, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(          0, spaceSize.y, spaceSize.z));
	geometry.vertices.push(new THREE.Vector3(          0, spaceSize.y,           0));

	rightSide = new THREE.Line(geometry, material);

	scene.add(rightSide);

	// add a point light
	pointLight = new THREE.PointLight(0xffffff);

	// set its position
	pointLight.position.x = x;
	pointLight.position.y = y;
	pointLight.position.z = z;

	// add to the scene
	scene.add(pointLight);

	sphereMaterial = new THREE.MeshLambertMaterial({ color: 0x0000cc });

}

// update particles position
function updateParticles() {

	end = xhrStream.responseText.indexOf(";", start);

	try {

		if (end > 0 && start >= 0) {
			particles.length = 0;
			particles = JSON.parse(xhrStream.responseText.substring(start, end));
		}

	} catch (e) {
		console.log(e);
	}

	start = end;
	start++;

	for (var i = 0, len = particles.length; i < len; i++) {

		if ( ! spheresMap.hasOwnProperty(particles[i].id)) {

			spheresMap[particles[i].id] = new THREE.Mesh(new THREE.SphereGeometry(particles[i].radius, 16, 16),sphereMaterial);

			spheresMap[particles[i].id].geometry.dynamic = true;
			spheresMap[particles[i].id].geometry.verticesNeedUpdate = true;
			spheresMap[particles[i].id].geometry.normalsNeedUpdate = true;

			scene.add(spheresMap[particles[i].id]);

		}

		spheresMap[particles[i].id].position.x = particles[i].pos.x;
		spheresMap[particles[i].id].position.y = particles[i].pos.y;
		spheresMap[particles[i].id].position.z = particles[i].pos.z;

	}

}

function animate() {

	requestAnimationFrame(animate);

	render();

}

function render() {

	camera.lookAt({ x: spaceSize.x/2, y: spaceSize.y/2, z: spaceSize.z/2 });

	updateParticles();

	renderer.render(scene, camera);

}

// Event handlers ---------------------------------------------------------- //

var isMouseDown = false;

var prevClickPos = { x: 0,  y: 0 };

// Mouse handler
function onMouseDown(event) {

	prevClickPos.x = event.x;
	prevClickPos.y = event.y;

	isMouseDown = true;

}

function onMouseMove(event) {

	var dx, dy;

	if (isMouseDown) {

		dx = event.x - prevClickPos.x;
		dy = event.y - prevClickPos.y;

		theta -= 0.005*dx;
		phi -= 0.005*dy;

		x = spaceSize.x/2 + r*Math.sin(theta)*Math.cos(phi);
		y = spaceSize.y/2 + r*Math.sin(theta)*Math.sin(phi);
		z = spaceSize.z/2 + r*Math.cos(theta);

		camera.position.x = x;
		camera.position.y = y;
		camera.position.z = z;

		pointLight.position.x = x;
		pointLight.position.y = y;
		pointLight.position.z = z;

		prevClickPos.x = event.x;
		prevClickPos.y = event.y;

	}

}

function onMouseUp(event) {

	isMouseDown = false;

}

function onKeyPress(event) {

	console.log(event);

}

function requestSettings() {

	var xhr = new XMLHttpRequest();

	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4 && xhr.status == 200) {
			onSettingsSuccess(xhr.responseText);
		}
	}

	xhr.open("GET", "/settings", true);
	xhr.send(null);
}

function onSettingsSuccess(responseText) {

	settings = JSON.parse(responseText);
	spaceSize = settings.simSpaceSize;

	initEnvironment();
	requestStream();
	animate();

}

function requestStream() {

	xhrStream = new XMLHttpRequest();

	xhrStream.onreadystatechange = function() {
		//onStreamLoad(xhrStream);
	}

	xhrStream.open("GET", "/simstream?rate=30", true);
	xhrStream.send(null);

}

function onStreamLoad(xhr) {

}

function resizeHandler(e) {

	view.width = document.getElementById("view").clientWidth;
	view.height = document.getElementById("view").clientHeight;

	var heights = { 
		top: document.getElementById("top").clientHeight, 
		bottom: document.getElementById("bottom").clientHeight
	};

	var sideWidth = document.getElementById("side").clientWidth - 1;

	document.getElementById("main").style.height = (window.innerHeight - heights.top - heights.bottom - 2)+"px";
	document.getElementById("view").style.width = (window.innerWidth - sideWidth - 1 - 20)+"px";

	renderer.setSize(view.width, view.height
		/* window.innerWidth, window.innerHeight */
	);

}

function loadHandler(e) {

	var heights = { 
		top : document.getElementById("top").clientHeight, 
		bottom : document.getElementById("bottom").clientHeight
	};

	var sideWidth = document.getElementById("side").clientWidth - 1;

	document.getElementById("main").style.height = (window.innerHeight - heights.top - heights.bottom - 2)+"px";
	document.getElementById("view").style.width = (window.innerWidth - sideWidth - 1 - 20)+"px";

	window.addEventListener("resize", resizeHandler, false);
	window.addEventListener("keypress", onKeyPress, false);

	document.getElementById("view").addEventListener("mousedown", onMouseDown, false);
	document.getElementById("view").addEventListener("mousemove", onMouseMove, false);
	document.getElementById("view").addEventListener("mouseup", onMouseUp, false);

	requestSettings();

}

window.addEventListener("load", loadHandler, false);