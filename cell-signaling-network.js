// Node.js script to populate a simulation domain with particles
// without overlapping
var fs = require('fs');

// Number of space cells per side
var Nx, Ny, Nz, N;

var simtime = 1000;

// Domain
var spaceSize = {
	"x": 1000,
	"y": 1000,
	"z": 1000
};

var spaceCellSize = 200;

var particlesLists = [];

var numParticles = 10000;

var particleRadius = 15;

var particleMass = 1;

var result = [];

var emitter = {
	"x": 500,
	"y": 250,
	"z": 500,
	"r": 100,
	"m": 1000
};

var receiver = {
	"x": 500,
	"y": 750,
	"z": 500,
	"r": 100,
	"m": 1000
};

var e = 0.0001;

function init() {

	Nx = Math.ceil(spaceSize.x/spaceCellSize);
	Ny = Math.ceil(spaceSize.y/spaceCellSize);
	Nz = Math.ceil(spaceSize.z/spaceCellSize);

	N = Nx*Ny*Nz;

	for (i = 0; i < N; i++) {
		particlesLists[i] = [];
	}

	particlesLists[getSpaceCellIndex(emitter)].push(emitter);
	particlesLists[getSpaceCellIndex(receiver)].push(receiver);

}

function run() {

	var count = 0;

	var index = 0;

	var i, j, k;

	var particleCellList;

	var particles = [];

	var overlap = false;

	while (count < numParticles) {

		particles.length = 0;
		overlap = false;

		var particle = {
			"x": particleRadius + e + (spaceSize.x-2*particleRadius)*Math.random(),
			"y": particleRadius + e + (spaceSize.y-2*particleRadius)*Math.random(),
			"z": particleRadius + e + (spaceSize.z-2*particleRadius)*Math.random(),
			"r": particleRadius
		};

		i = Math.floor(particle.x/spaceSize.x);
		j = Math.floor(particle.y/spaceSize.y);
		k = Math.floor(particle.z/spaceSize.z);

		for (var a = -1; a <= 1; a++) {
		for (var b = -1; b <= 1; b++) {
		for (var c = -1; c <= 1; c++) {

			if (belongsToSpace(i+a, j+b, k+c)) {
				particles = particles.concat(particlesLists[i*Ny*Nz + j*Nz + k]);
			}

		}}}

		for (var p = 0; p < particles.length; p++) {

			overlap = checkOverlap(particle, particles[p]);

			if (overlap) {
				break;
			}

		}

		if ( ! overlap) {
			particlesLists[i*Ny*Nz + j*Nz + k].push(particle);
			result.push(particle);
			console.log(
				"particle["+count+"]:" +
				"  x = " + particle.x +
				", y = " + particle.y + 
				", z = " + particle.z);
			count++;
		}
	}

	console.log("Completed!");
}

function belongsToSpace(i, j, k) {
	return ((0 <= i && i < Nx) &&
		(0 <= j && j < Ny) &&
		(0 <= k && k < Nz));
}

function checkOverlap(pa, pb) {

	var dx = (pa.x-pb.x);
	var dy = (pa.y-pb.y);
	var dz = (pa.z-pb.z);

	return  (Math.sqrt(dx*dx + dy*dy + dz*dz) < pa.r + pb.r ? true : false);
}

function getSpaceCellIndex(p) {

	var i = Math.floor(p.x/spaceSize.x), 
		j = Math.floor(p.y/spaceSize.y), 
		k = Math.floor(p.z/spaceSize.z);

	return i*Ny*Nz + j*Nz + k; 
}

function writeIniFile() {

	fs.open("networks/cell-signaling.ini", 'w', function(err, fd) {

		fs.writeSync(fd, "[Config cell-signaling]\n");
		fs.writeSync(fd, "network = domain\n");
		fs.writeSync(fd, "sim-time-limit = "+simtime+"s\n");
		fs.writeSync(fd, "\n");

		fs.writeSync(fd, "domain.spaceSizeX = "+spaceSize.x+"\n");
		fs.writeSync(fd, "domain.spaceSizeY = "+spaceSize.y+"\n");
		fs.writeSync(fd, "domain.spaceSizeZ = "+spaceSize.z+"\n");
		fs.writeSync(fd, "domain.numberOfInitialMolecules = "+numParticles+"\n");
		fs.writeSync(fd, "\n");

		fs.writeSync(fd, "domain.manager.enableWebServer = 0\n");
		fs.writeSync(fd, "domain.manager.tkRefreshRate = 100\n");
		fs.writeSync(fd, "domain.manager.statsRefreshRate = 1000\n");
		fs.writeSync(fd, "domain.manager.mode = 1\n");
		fs.writeSync(fd, "\n");

		fs.writeSync(fd, "domain.cell[0].mobility.radius = "+emitter.r+"\n");
		fs.writeSync(fd, "domain.cell[0].mobility.mass = "+emitter.m+"\n");
		fs.writeSync(fd, "domain.cell[0].mobility.xpos = "+emitter.x+"\n");
		fs.writeSync(fd, "domain.cell[0].mobility.ypos = "+emitter.y+"\n");
		fs.writeSync(fd, "domain.cell[0].mobility.zpos = "+emitter.z+"\n");
		fs.writeSync(fd, "domain.cell[0].mobility.vx = 0\n");
		fs.writeSync(fd, "domain.cell[0].mobility.vy = 0\n");
		fs.writeSync(fd, "domain.cell[0].mobility.vz = 0\n");
		fs.writeSync(fd, "domain.cell[0].emitter.enabled = true\n");
		fs.writeSync(fd, "domain.cell[0].receiver.enabled = false\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionStart = 10\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionDuration = 200\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionRate = 500\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionParticleRadius = 0.2\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionParticleMass = 1\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionTimeToLive = "+simtime+"\n");
		fs.writeSync(fd, "domain.cell[0].emitter.emissionBoundariesMode = 2\n");
		fs.writeSync(fd, "\n");

		fs.writeSync(fd, "domain.cell[1].mobility.radius = "+receiver.r+"\n");
		fs.writeSync(fd, "domain.cell[1].mobility.mass = "+receiver.m+"\n");
		fs.writeSync(fd, "domain.cell[1].mobility.xpos = "+receiver.x+"\n");
		fs.writeSync(fd, "domain.cell[1].mobility.ypos = "+receiver.y+"\n");
		fs.writeSync(fd, "domain.cell[1].mobility.zpos = "+receiver.z+"\n");
		fs.writeSync(fd, "domain.cell[1].mobility.vx = 0\n");
		fs.writeSync(fd, "domain.cell[1].mobility.vy = 0\n");
		fs.writeSync(fd, "domain.cell[1].mobility.vz = 0\n");
		fs.writeSync(fd, "domain.cell[1].receiver.enabled = true\n");
		fs.writeSync(fd, "domain.cell[1].emitter.enabled = false\n");
		fs.writeSync(fd, "domain.cell[1].receiver.statsRefreshRate = 1000\n");
		fs.writeSync(fd, "\n");

		fs.writeSync(fd, "domain.molecule[*].mass = "+particleMass+"\n");
		fs.writeSync(fd, "domain.molecule[*].radius = "+particleRadius+"\n");
		fs.writeSync(fd, "\n");

		fs.writeSync(fd, "domain.molecule[*].vx = normal(0,1)\n");
		fs.writeSync(fd, "domain.molecule[*].vy = normal(0,1)\n");
		fs.writeSync(fd, "domain.molecule[*].vz = normal(0,1)\n");
		fs.writeSync(fd, "\n");

		for (var i = 0; i < result.length; i++) {
			fs.writeSync(fd, "domain.molecule["+i+"].xpos = "+result[i].x+"\n");
			fs.writeSync(fd, "domain.molecule["+i+"].ypos = "+result[i].y+"\n");
			fs.writeSync(fd, "domain.molecule["+i+"].zpos = "+result[i].z+"\n");
			fs.writeSync(fd, "\n");
		}

	});
}

init();

run();

writeIniFile();