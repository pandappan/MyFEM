# MyFEM JSON Input Schema

## Top-level

| Field | Type | Required | Description |
|---|---|---|---|
| title | string | yes | Model description |
| dimension | int (2 or 3) | yes | Spatial dimension |
| nodes | array | yes | Node list |
| materials | array | yes | Material list |
| element_groups | array | yes | Element group list |
| boundary_conditions | object | yes | BC definitions |
| loads | object | yes | Load definitions |

## nodes[]

| Field | Type | Required | Default | Description |
|---|---|---|---|---|
| id | int | yes | — | Unique node id (1-based) |
| x | double | yes | — | X coordinate |
| y | double | yes | — | Y coordinate |
| z | double | no | 0.0 | Z coordinate (omit for 2D) |

## materials[]

| Field | Type | Required | Default | Description |
|---|---|---|---|---|
| id | int | yes | — | Unique material id |
| type | string | yes | — | "elastic_isotropic" |
| E | double | yes | — | Young's modulus |
| nu | double | yes | — | Poisson's ratio |
| rho | double | no | 0.0 | Density |
| thickness | double | no | 1.0 | For 2D plane elements |
| area | double | no | 1.0 | For bar elements |

## element_groups[]

| Field | Type | Required | Description |
|---|---|---|---|
| type | string | yes | "Bar3D" / "Q4_PS" / "Q4_PE" / "H8" |
| material | int | yes | Material id reference |
| elements | array | yes | Element list |

### element_groups[].elements[]

| Field | Type | Required | Description |
|---|---|---|---|
| id | int | yes | Unique element id |
| nodes | int[] | yes | Node id list (element-type-specific order) |

## boundary_conditions

### .fixed[]
| Field | Type | Required | Description |
|---|---|---|---|
| nodes | int[] | yes | Node ids to constrain |
| dof | string[] | yes | Subset of ["x","y","z"] |

### .prescribed[]
| Field | Type | Required | Description |
|---|---|---|---|
| node | int | yes | Node id |
| dof | string | yes | "x"/"y"/"z" |
| value | double | yes | Prescribed displacement |

## loads

### .concentrated[]
| Field | Type | Required | Description |
|---|---|---|---|
| node | int | yes | Node id |
| dof | string | yes | "x"/"y"/"z" |
| value | double | yes | Force value (per node) |

### .surface[]
| Field | Type | Required | Description |
|---|---|---|---|
| element | int | yes | Element id |
| face | int | yes | Face id (1-based, see element face convention) |
| dof | string | yes | "x"/"y"/"z" |
| value | double | yes | Distributed load per unit length/area |

### .body_force
Array of 3 doubles: [gx, gy, gz]. Optional, default [0,0,0].