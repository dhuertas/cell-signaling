// Simulation environment -------------------------------------------------- //
var settings;

var spaceSize;

var particles = [];

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

// XHR stream object and variables to parse the received data so far
var stream;
var rate;
var start, end;

var view = {
	width: 0,
	height: 0
};

// Stats
var stats = {
	simulationTime: 0,
	numberOfParticles: 0,
	spaceSizeX: 0,
	spaceSizeY: 0,
	spaceSizeZ:0
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

	var tmp;

	tmp = stream.responseText.indexOf(";", start);

	if (tmp == -1) {
		// We have reached the end of the received data so far
		// console.log("end of stream reached");
		return;
	}

	end = tmp;

	try {

		if ( 0 <= start && start < end && end > 0) {

			particles.length = 0;
			particles = JSON.parse(stream.responseText.substring(start, end));

			stats.numberOfParticles = particles.length;

			start = end;
			start++;

		}

		for (var i = 0, len = particles.length; i < len; i++) {

			if (particles[i].hasOwnProperty("stats")) {
				// Received a stats JSON object update
				stats.simulationTime = particles[i].st;
				continue;
			}

			if ( ! spheresMap.hasOwnProperty(particles[i].id)) {

				spheresMap[particles[i].id] = new THREE.Mesh(
					new THREE.SphereGeometry(particles[i].radius, 16, 16),
					sphereMaterial
				);

				spheresMap[particles[i].id].geometry.dynamic = true;
				spheresMap[particles[i].id].geometry.verticesNeedUpdate = true;
				spheresMap[particles[i].id].geometry.normalsNeedUpdate = true;

				scene.add(spheresMap[particles[i].id]);

			}

			spheresMap[particles[i].id].position.x = particles[i].pos.x;
			spheresMap[particles[i].id].position.y = particles[i].pos.y;
			spheresMap[particles[i].id].position.z = particles[i].pos.z;
			spheresMap[particles[i].id].updated = true;

		}

		// Remove the spheres that have not been updated
		for (var sphere in spheresMap) {
			if (spheresMap[sphere].updated == false) {
				// This sphere has not been updated, remove it
				delete spheresMap[sphere];
			} else {
				spheresMap[sphere].updated = false;
			}
		}

	} catch (e) {
		console.log(e);
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

	// Update stats every 2 secs
	setInterval(updateStatistics, 2000);

	initEnvironment();
	requestStream();
	animate();

}

function requestStream() {

	rate = 30;
	start = 0;
	end = 0;

	stream = new XMLHttpRequest();

	stream.onreadystatechange = function() {
		//onStreamLoad(stream);
	}

	stream.open("GET", "/simstream?rate=" + rate, true);
	stream.send(null);

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

	if (typeof renderer !== "undefined") {
		renderer.setSize(view.width, view.height
			/* window.innerWidth, window.innerHeight */
		);
	}

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

	document.getElementById("start").addEventListener("click", start, false);
	document.getElementById("stop").addEventListener("click", stop, false);

	requestSettings();

}

// Statistics -------------------------------------------------------------- //

function updateStatistics() {

	// Get the DOM elements and update them
	document.getElementById("number-of-particles").innerText = stats.numberOfParticles;
	document.getElementById("sim-time").innerText = stats.simulationTime;
	document.getElementById("space-size-x").innerText = spaceSize.x;
	document.getElementById("space-size-y").innerText = spaceSize.y;
	document.getElementById("space-size-z").innerText = spaceSize.z;

}

// Actions ----------------------------------------------------------------- //

function start() {

	var xhr = new XMLHttpRequest();

	xhr.onreadystatechange = function() {}

	xhr.open("GET", "action?cmd=start", true);
	xhr.send(null);

}

function stop() {

	var xhr = new XMLHttpRequest();

	xhr.onreadystatechange = function() {}

	xhr.open("GET", "action?cmd=stop", true);
	xhr.send(null);

}

window.addEventListener("load", loadHandler, false);