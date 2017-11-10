! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

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
