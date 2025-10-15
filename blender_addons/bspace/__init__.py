#####################################################################################################################
# Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#####################################################################################################################


bl_info = {
    "name": "BSpace",
    "author": "Milan Jaros, Petr Strakos, Lubomir Riha",
    "description": "",
    "blender": (4, 2, 0),
    "version": (0, 3, 0),
    "location": "View3D > Sidebar > BSpace Tab",
    "warning": "",
    "category": "3D View"
}
#####################################################################################################################

def register():
    from . import bspace_pref
    from . import bspace_panel

    bspace_pref.register()
    bspace_panel.register()

def unregister():
    from . import bspace_pref
    from . import bspace_panel
    
    try:        
        bspace_pref.unregister()
        bspace_panel.unregister()

    except RuntimeError:
        pass 