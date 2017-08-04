
!> Linked list interface block

!> Node of a linked list
type :: node_t
  integer           :: key
  type(LISTED_TYPE) :: element

  type(node_t), pointer  :: next => NULL()
end type

!> Registry type
type :: registry_t
  logical               :: l_init = .false.
  integer               :: count  = 0
  type(node_t), pointer :: head   => NULL()

contains
  procedure :: init => init_
  procedure :: finalize => finalize_
  procedure :: add => add_
  procedure :: get => get_
  procedure :: remove => remove_
end type
