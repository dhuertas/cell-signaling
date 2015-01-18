//  This file is part of the Cell-Signaling project. Cell-Signaling is an
//  Omnet++ project to simulate cell signaling communications.
//  Copyright (C) 2014  Daniel Huertas
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package cellsignaling.networks;

import cellsignaling.src.Cell;
import cellsignaling.src.Manager;
import cellsignaling.src.Molecule;

network domain
{
  @display("bgb=400,300");

  double spaceSizeX = default(1000);
  double spaceSizeY = default(1000);
  double spaceSizeZ = default(1000);

  int numberOfInitialMolecules = default(0);

  submodules:
    manager: Manager {
      @display("p=-30,-30;i=block/cogwheel");
    }

    cell[2]: Cell;
    molecule[numberOfInitialMolecules]: Molecule;

  connections allowunconnected:

}
